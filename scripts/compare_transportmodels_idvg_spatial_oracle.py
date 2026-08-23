#!/usr/bin/env python3
"""Compare five-bias Sentaurus/Vela TransportModels DG spatial profiles."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
SENTAURUS_MANIFEST = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "idvg_spatial_oracle_20260821/spatial_oracle_manifest.json"
)
VELA_MANIFEST = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "idvg_spatial_oracle_2026-08-21/vela/spatial_oracle_manifest.json"
)
OUTPUT_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "idvg_spatial_oracle_2026-08-21/comparison"
)
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_idvg_spatial_oracle_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_idvg_spatial_oracle_2026-08-21.md"
PROFILE_TARGETS = (
    ("source_end", -0.020),
    ("channel_mid", 0.000),
    ("drain_end", 0.020),
)
MAX_DEPTH_UM = 0.020


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = int(math.floor(position)), int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def read_scalar_csv(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["component0"]) for row in csv.DictReader(handle)}


def read_vector_csv(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
            for row in csv.DictReader(handle)
        }


def read_nodes(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
            for row in csv.DictReader(handle)
        }


def parse_vtk(path: Path) -> tuple[list[tuple[float, float]], dict[str, list[Any]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points_index = next(i for i, line in enumerate(lines) if line.startswith("POINTS "))
    point_count = int(lines[points_index].split()[1])
    points = []
    for line in lines[points_index + 1 : points_index + 1 + point_count]:
        values = line.split()
        points.append((float(values[0]), float(values[1])))
    point_data_index = next(i for i, line in enumerate(lines) if line.startswith("POINT_DATA "))
    data_count = int(lines[point_data_index].split()[1])
    fields: dict[str, list[Any]] = {}
    index = point_data_index + 1
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
            continue
        if parts[0] == "SCALARS":
            name = parts[1]
            index += 2  # Skip SCALARS and LOOKUP_TABLE.
            fields[name] = [float(lines[index + offset].split()[0]) for offset in range(data_count)]
            index += data_count
        elif parts[0] == "VECTORS":
            name = parts[1]
            index += 1
            fields[name] = [
                tuple(float(value) for value in lines[index + offset].split()[:3])
                for offset in range(data_count)
            ]
            index += data_count
        else:
            index += 1
    return points, fields


def read_vela_state(path: Path) -> dict[int, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): {key: float(value) for key, value in row.items() if key != "node_id"} for row in csv.DictReader(handle)}


def sentaurus_fields(root: Path) -> dict[str, dict[int, Any]]:
    fields = root / "fields"
    return {
        "qn_V": read_scalar_csv(fields / "eQuantumPotential_region3.csv"),
        "n_cm3": read_scalar_csv(fields / "eDensity_region3.csv"),
        "phin_V": read_scalar_csv(fields / "eQuasiFermiPotential_region3.csv"),
        "enormal_V_cm": read_scalar_csv(fields / "eEnormal_region3.csv"),
        "mobility_cm2_V_s": read_scalar_csv(fields / "eMobility_region3.csv"),
        "grad_phin_V_cm": read_vector_csv(fields / "eGradQuasiFermi_region3.csv"),
    }


def profile_nodes(nodes: dict[int, tuple[float, float]], available: set[int]) -> dict[str, dict[str, Any]]:
    silicon_y = {y for node, (x, y) in nodes.items() if node in available and x >= -1.0e-12}
    profiles: dict[str, dict[str, Any]] = {}
    for name, target in PROFILE_TARGETS:
        selected_y = min(silicon_y, key=lambda value: abs(value - target))
        ids = sorted(
            (
                node
                for node, (x, y) in nodes.items()
                if node in available
                and abs(y - selected_y) < 1.0e-12
                and -1.0e-12 <= x <= MAX_DEPTH_UM + 1.0e-12
            ),
            key=lambda node: nodes[node][0],
        )
        if len(ids) < 10:
            raise RuntimeError(f"Profile {name} has only {len(ids)} nodes")
        profiles[name] = {"target_y_um": target, "selected_y_um": selected_y, "node_ids": ids}
    return profiles


def field_metrics(rows: list[dict[str, Any]], error_key: str) -> dict[str, float]:
    values = [float(row[error_key]) for row in rows]
    return {
        "median": percentile(values, 0.5),
        "p95": percentile(values, 0.95),
        "max": max(values),
    }


def make_plots(long_rows: list[dict[str, Any]], summary: list[dict[str, Any]]) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    chosen_bias = 1.0
    selected = [row for row in long_rows if abs(row["gate_bias_V"] - chosen_bias) < 1.0e-12]
    fields = (
        ("qn", "Quantum potential Qn (V)", False),
        ("n", "Electron density n (cm⁻³)", True),
        ("phin", "Electron quasi-Fermi φn (V)", False),
        ("enormal", "|Enormal| (V/cm)", True),
        ("mobility", "Electron mobility (cm²/Vs)", False),
        ("grad_phin", "|∇φn| (V/cm)", True),
    )
    colors = {"source_end": "tab:blue", "channel_mid": "tab:green", "drain_end": "tab:red"}
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.3))
    for axis, (key, ylabel, logarithmic) in zip(axes.flat, fields):
        for profile, color in colors.items():
            rows = [row for row in selected if row["profile"] == profile]
            axis.plot([row["depth_nm"] for row in rows], [row[f"sentaurus_{key}"] for row in rows], "-", color=color, label=f"Sentaurus {profile}")
            axis.plot([row["depth_nm"] for row in rows], [row[f"vela_{key}"] for row in rows], "--", color=color, label=f"Vela {profile}")
        if logarithmic:
            axis.set_yscale("log")
        axis.set_xlabel("Depth into silicon (nm)")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("TransportModels spatial profiles at Vg=1.00 V, Vd=1.10 V")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    profile_png = OUTPUT_ROOT / "idvg_spatial_profiles_vg1p00.png"
    profile_svg = OUTPUT_ROOT / "idvg_spatial_profiles_vg1p00.svg"
    fig.savefig(profile_png, dpi=180)
    fig.savefig(profile_svg)
    plt.close(fig)

    metric_specs = (
        ("qn_mV", "Qn p95 (mV)"),
        ("n_dex", "n p95 (dex)"),
        ("phin_mV", "φn p95 (mV)"),
        ("enormal_dex", "Enormal p95 (dex)"),
        ("mobility_dex", "mobility p95 (dex)"),
        ("grad_phin_dex", "|∇φn| p95 (dex)"),
    )
    biases = sorted({float(row["gate_bias_V"]) for row in summary})
    values = np.array(
        [
            [next(row[f"{key}_p95"] for row in summary if abs(row["gate_bias_V"] - bias) < 1.0e-12) for bias in biases]
            for key, _ in metric_specs
        ]
    )
    normalized = values / np.maximum(values.max(axis=1, keepdims=True), 1.0e-30)
    fig, axis = plt.subplots(figsize=(9.0, 5.2))
    image = axis.imshow(normalized, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=1.0)
    axis.set_xticks(range(len(biases)), [f"{bias:.2f}" for bias in biases])
    axis.set_yticks(range(len(metric_specs)), [label for _, label in metric_specs])
    axis.set_xlabel("Gate voltage Vg (V)")
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            axis.text(column_index, row_index, f"{values[row_index, column_index]:.3g}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=axis, label="Row-normalized error")
    axis.set_title("Five-bias spatial p95 errors across three normal profiles")
    fig.tight_layout()
    heatmap_png = OUTPUT_ROOT / "idvg_spatial_error_heatmap.png"
    heatmap_svg = OUTPUT_ROOT / "idvg_spatial_error_heatmap.svg"
    fig.savefig(heatmap_png, dpi=180)
    fig.savefig(heatmap_svg)
    plt.close(fig)
    return {
        "profile_png": str(profile_png.resolve()),
        "profile_svg": str(profile_svg.resolve()),
        "heatmap_png": str(heatmap_png.resolve()),
        "heatmap_svg": str(heatmap_svg.resolve()),
    }


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels five-bias spatial oracle comparison",
        "",
        f"Status: **{report['status']}**. Node mapping: **{report['node_mapping_status']}**.",
        "",
        "Profiles follow the silicon normal from the Si/SiO2 interface to 20 nm depth at the source end, channel midpoint, and drain end.",
        "",
        "| Vg (V) | Qn p95 (mV) | n p95 (dex) | φn p95 (mV) | Enormal p95 (dex) | μn p95 (dex) | |∇φn| p95 (dex) |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in report["bias_summary"]:
        lines.append(
            f"| {row['gate_bias_V']:.2f} | {row['qn_mV_p95']:.6g} | {row['n_dex_p95']:.6g} | "
            f"{row['phin_mV_p95']:.6g} | {row['enormal_dex_p95']:.6g} | "
            f"{row['mobility_dex_p95']:.6g} | {row['grad_phin_dex_p95']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Profile localization",
            "",
            "| Vg (V) | Worst Qn profile / p95 (mV) | Worst n profile / p95 (dex) | Worst φn profile / p95 (mV) |",
            "|---:|---:|---:|---:|",
        ]
    )
    for bias in sorted({row["gate_bias_V"] for row in report["profile_summary"]}):
        rows = [row for row in report["profile_summary"] if row["gate_bias_V"] == bias]
        qn = max(rows, key=lambda row: row["qn_mV_p95"])
        density = max(rows, key=lambda row: row["n_dex_p95"])
        phin = max(rows, key=lambda row: row["phin_mV_p95"])
        lines.append(
            f"| {bias:.2f} | {qn['profile']} / {qn['qn_mV_p95']:.6g} | "
            f"{density['profile']} / {density['n_dex_p95']:.6g} | "
            f"{phin['profile']} / {phin['phin_mV_p95']:.6g} |"
        )
    lines.extend(
        [
            "",
            f"Qn criterion (`p95 <= 20 mV`): **{'pass' if report['acceptance']['qn_p95_le_20mV'] else 'fail'}**.",
            f"Electron-density criterion (`p95 <= 0.2 dex`): **{'pass' if report['acceptance']['n_p95_le_0p2dex'] else 'fail'}**.",
            "",
            f"Profile figure: `{report['artifacts']['profile_png']}`",
            f"Error heatmap: `{report['artifacts']['heatmap_png']}`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        assert sha256(Path(report["artifacts"]["long_csv"])) == report["artifacts"]["long_csv_sha256"]
        assert report["node_mapping_status"] == "exact"
        print("TransportModels Id-Vg spatial oracle check: PASS")
        return 0

    sent_manifest = json.loads(SENTAURUS_MANIFEST.read_text(encoding="utf-8"))
    vela_manifest = json.loads(VELA_MANIFEST.read_text(encoding="utf-8"))
    sent_states = {round(float(row["gate_bias_V"]), 10): row for row in sent_manifest["states"]}
    vela_states = {round(float(row["gate_bias_V"]), 10): row for row in vela_manifest["states"]}
    if sent_states.keys() != vela_states.keys():
        raise RuntimeError("Sentaurus/Vela bias sets differ")

    long_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    profile_summaries: list[dict[str, Any]] = []
    profile_metadata = None
    exact_mapping = True
    for bias in sorted(sent_states):
        sent_root = Path(sent_states[bias]["export_dir"])
        nodes = read_nodes(sent_root / "nodes.csv")
        sent = sentaurus_fields(sent_root)
        vtk_points, vtk = parse_vtk(Path(vela_states[bias]["vtk"]))
        vela_state = read_vela_state(Path(vela_states[bias]["state_csv"]))
        if len(vtk_points) != len(nodes):
            exact_mapping = False
        else:
            for node, point in enumerate(vtk_points):
                if node not in nodes or max(abs(point[0] - nodes[node][0]), abs(point[1] - nodes[node][1])) > 1.0e-12:
                    exact_mapping = False
                    break
        available = set(sent["n_cm3"]) & set(vela_state)
        profiles = profile_nodes(nodes, available)
        if profile_metadata is None:
            profile_metadata = profiles
        bias_rows: list[dict[str, Any]] = []
        for profile_name, profile in profiles.items():
            for node in profile["node_ids"]:
                sent_grad = math.hypot(*sent["grad_phin_V_cm"][node])
                vela_enormal = abs(float(vtk["ElectricFieldVector"][node][0]))
                sent_enormal = abs(sent["enormal_V_cm"][node])
                # The unit-scaled VTK writer exports mobility in cm^2/(V*s),
                # matching the Sentaurus eMobility dataset directly.
                vela_mobility = float(vtk["ElectronMobility"][node])
                sent_mobility = abs(sent["mobility_cm2_V_s"][node])
                vela_grad = abs(float(vtk["ElectronHighFieldDrive"][node]))
                sent_density = max(sent["n_cm3"][node], 1.0)
                vela_density = max(vela_state[node]["electrons_m3"] / 1.0e6, 1.0)
                row = {
                    "gate_bias_V": bias,
                    "profile": profile_name,
                    "node_id": node,
                    "depth_nm": max(nodes[node][0], 0.0) * 1.0e3,
                    "sentaurus_qn": sent["qn_V"][node],
                    "vela_qn": vela_state[node]["electron_quantum_potential_V"],
                    "sentaurus_n": sent_density,
                    "vela_n": vela_density,
                    "sentaurus_phin": sent["phin_V"][node],
                    "vela_phin": float(vtk["ElectronQuasiFermi"][node]),
                    "sentaurus_enormal": sent_enormal,
                    "vela_enormal": vela_enormal,
                    "sentaurus_mobility": sent_mobility,
                    "vela_mobility": vela_mobility,
                    "sentaurus_grad_phin": sent_grad,
                    "vela_grad_phin": vela_grad,
                }
                row.update(
                    {
                        "qn_error_mV": abs(row["vela_qn"] - row["sentaurus_qn"]) * 1.0e3,
                        "n_error_dex": abs(math.log10(row["vela_n"]) - math.log10(row["sentaurus_n"])),
                        "phin_error_mV": abs(row["vela_phin"] - row["sentaurus_phin"]) * 1.0e3,
                        "enormal_error_dex": abs(math.log10(max(row["vela_enormal"], 1.0)) - math.log10(max(row["sentaurus_enormal"], 1.0))),
                        "mobility_error_dex": abs(math.log10(max(row["vela_mobility"], 1.0e-12)) - math.log10(max(row["sentaurus_mobility"], 1.0e-12))),
                        "grad_phin_error_dex": abs(math.log10(max(row["vela_grad_phin"], 1.0)) - math.log10(max(row["sentaurus_grad_phin"], 1.0))),
                    }
                )
                bias_rows.append(row)
                long_rows.append(row)
        for profile_name in profiles:
            selected_rows = [row for row in bias_rows if row["profile"] == profile_name]
            profile_summaries.append(
                {
                    "gate_bias_V": bias,
                    "profile": profile_name,
                    **{f"qn_mV_{key}": value for key, value in field_metrics(selected_rows, "qn_error_mV").items()},
                    **{f"n_dex_{key}": value for key, value in field_metrics(selected_rows, "n_error_dex").items()},
                    **{f"phin_mV_{key}": value for key, value in field_metrics(selected_rows, "phin_error_mV").items()},
                    **{f"enormal_dex_{key}": value for key, value in field_metrics(selected_rows, "enormal_error_dex").items()},
                    **{f"mobility_dex_{key}": value for key, value in field_metrics(selected_rows, "mobility_error_dex").items()},
                    **{f"grad_phin_dex_{key}": value for key, value in field_metrics(selected_rows, "grad_phin_error_dex").items()},
                }
            )
        summaries.append(
            {
                "gate_bias_V": bias,
                **{f"qn_mV_{key}": value for key, value in field_metrics(bias_rows, "qn_error_mV").items()},
                **{f"n_dex_{key}": value for key, value in field_metrics(bias_rows, "n_error_dex").items()},
                **{f"phin_mV_{key}": value for key, value in field_metrics(bias_rows, "phin_error_mV").items()},
                **{f"enormal_dex_{key}": value for key, value in field_metrics(bias_rows, "enormal_error_dex").items()},
                **{f"mobility_dex_{key}": value for key, value in field_metrics(bias_rows, "mobility_error_dex").items()},
                **{f"grad_phin_dex_{key}": value for key, value in field_metrics(bias_rows, "grad_phin_error_dex").items()},
            }
        )

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    long_csv = OUTPUT_ROOT / "idvg_spatial_profiles_long.csv"
    with long_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(long_rows[0]))
        writer.writeheader()
        writer.writerows(long_rows)
    artifacts = make_plots(long_rows, summaries)
    artifacts.update({"long_csv": str(long_csv.resolve()), "long_csv_sha256": sha256(long_csv)})
    acceptance = {
        "qn_p95_le_20mV": max(row["qn_mV_p95"] for row in summaries) <= 20.0,
        "n_p95_le_0p2dex": max(row["n_dex_p95"] for row in summaries) <= 0.2,
        "source_and_mid_qn_p95_le_20mV": max(
            row["qn_mV_p95"] for row in profile_summaries if row["profile"] != "drain_end"
        ) <= 20.0,
        "source_and_mid_n_p95_le_0p2dex": max(
            row["n_dex_p95"] for row in profile_summaries if row["profile"] != "drain_end"
        ) <= 0.2,
    }
    report = {
        "schema": "vela.transportmodels.idvg_spatial_oracle.v1",
        "as_of": "2026-08-21",
        "status": "complete",
        "node_mapping_status": "exact" if exact_mapping else "mismatch",
        "fixed_drain_bias_V": 1.1,
        "profiles": profile_metadata,
        "bias_summary": summaries,
        "profile_summary": profile_summaries,
        "acceptance": acceptance,
        "sentaurus_manifest": str(SENTAURUS_MANIFEST.resolve()),
        "vela_manifest": str(VELA_MANIFEST.resolve()),
        "artifacts": artifacts,
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "node_mapping_status": report["node_mapping_status"], "acceptance": acceptance, "bias_summary": summaries}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
