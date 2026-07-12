#!/usr/bin/env python3
"""Run PN2D coarse -18 V quasi-Fermi cap/warm-start branch diagnostics."""

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
BASE_REPORT = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "reports"
    / "coarse_vm_vector_compare"
)
BASE_CONFIG = BASE_REPORT / "simulation_coarse_previous_full20_aligned.json"
SENTAURUS_EXPORT = BASE_REPORT / "sentaurus_multibias"
RUNNER = REPO / "build-release" / ("vela_example_runner.exe" if os.name == "nt" else "vela_example_runner")
OUT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "reports"
    / "qf_cap_warmstart_branch_20260705"
)

CASES = [
    ("A", "warm_start=true, cap=0.1", True, 0.1, True, 1.0),
    ("B", "warm_start=true, cap=0.0", True, 0.0, True, 1.0),
    ("C", "warm_start=false, cap=0.0", False, 0.0, True, 1.0),
    ("D", "warm_start=true, cap=0.0, line_search=false", True, 0.0, False, 1.0),
]

BIAS_LADDER = [0.0, -1.0, -5.0, -10.0, -16.0, -18.0]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return [
            {str(k).strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            for row in csv.DictReader(fh)
        ]


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fnum(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def median(values: list[float]) -> float | None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def close(a: float | None, b: float, tol: float = 1.0e-6) -> bool:
    return a is not None and abs(a - b) <= tol


def configure_case(case_dir: Path, tag: str, warm: bool, cap: float, line_search: bool, damping: float) -> Path:
    cfg = read_json(BASE_CONFIG)
    solver = cfg.setdefault("solver", {})
    sweep = cfg.setdefault("sweep", {})
    solver["warm_start"] = warm
    solver["quasi_fermi_update_limit_V"] = cap
    solver["line_search"] = line_search
    solver["damping_factor"] = damping
    handoff = solver.setdefault("handoff", {})
    handoff["fallback"] = "none"
    handoff["require_gummel_convergence"] = False
    handoff["gummel_max_iter"] = 0
    handoff["newton_max_iter"] = solver.get("max_iter", handoff.get("newton_max_iter", 40))

    cfg["output_csv"] = f"coarse_qf_case_{tag}.csv"
    sweep["bias_points"] = BIAS_LADDER
    sweep["start"] = 0.0
    sweep["stop"] = -18.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = f"vtk_case_{tag}/dc_sweep"
    sweep["write_state_file"] = f"coarse_qf_case_{tag}_last_state.csv"
    sweep["diagnostics"] = {
        "sg_avalanche_edges": {"enabled": True, "csv_file": str(case_dir / f"sg_avalanche_edges_case_{tag}.csv")},
        "continuity_balance": {"enabled": True, "contacts": ["Anode", "Cathode"], "csv_file": str(case_dir / f"continuity_balance_case_{tag}.csv")},
        "newton_history": {"enabled": True, "csv_file": str(case_dir / f"newton_history_case_{tag}.csv")},
    }
    config_path = case_dir / f"simulation_case_{tag}.json"
    write_json(config_path, cfg)
    return config_path


def run_case(case_dir: Path, config_path: Path) -> tuple[int, str, str]:
    proc = subprocess.run(
        [str(RUNNER), "--config", str(config_path)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    (case_dir / "runner_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (case_dir / "runner_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    return proc.returncode, proc.stdout, proc.stderr


def import_compare_module():
    sys.path.insert(0, str(REPO / "scripts"))
    import run_pn2d_coarse7x3_previous_full20_compare as compare

    return compare


def build_compare_outputs(case_dir: Path, tag: str) -> tuple[Path, Path]:
    compare = import_compare_module()
    sentaurus_exports = {-18.0: SENTAURUS_EXPORT / "sentaurus_-18v"}
    vtk_dir = case_dir / f"vtk_case_{tag}"
    vela_vtks = compare.discover_vela_vtks(vtk_dir)
    node_rows = compare.node_field_compare(
        sentaurus_exports=sentaurus_exports,
        vela_vtks=vela_vtks,
        biases=[-18.0],
        exact_vela_bias=True,
    )
    node_csv = case_dir / f"coarse_node_field_compare_case_{tag}.csv"
    compare.write_csv_rows(node_csv, node_rows, [
        "bias_V", "vela_bias_V", "node_id", "x_um", "y_um", "quantity", "sentaurus_field",
        "sentaurus_value", "vela_field", "vela_value_scaled_to_sentaurus_units",
        "comparison_basis", "vela_over_sentaurus", "diff", "abs_diff",
    ])
    support_rows = compare.current_support_compare(
        sg_edges_csv=case_dir / f"sg_avalanche_edges_case_{tag}.csv",
        sentaurus_exports=sentaurus_exports,
    )
    support_csv = case_dir / f"coarse_current_support_compare_case_{tag}.csv"
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
    return node_csv, support_csv


def branch_summary_from_node_csv(path: Path) -> dict[str, Any]:
    by_node: dict[int, dict[str, float]] = {}
    for row in read_rows(path):
        if not close(fnum(row.get("bias_V")), -18.0):
            continue
        node = row.get("node_id")
        quantity = row.get("quantity", "").strip()
        if not node:
            continue
        item = by_node.setdefault(int(node), {})
        sv = fnum(row.get("sentaurus_value"))
        vv = fnum(row.get("vela_value_scaled_to_sentaurus_units"))
        if quantity in {"potential", "electron_qf", "hole_qf", "electron_density", "hole_density"}:
            if sv is not None:
                item[f"sent_{quantity}"] = sv
            if vv is not None:
                item[f"vela_{quantity}"] = vv
    dpsi_phin: list[float] = []
    logn: list[float] = []
    for item in by_node.values():
        needed = ["sent_potential", "sent_electron_qf", "vela_potential", "vela_electron_qf"]
        if all(k in item for k in needed):
            dpsi_phin.append((item["vela_potential"] - item["vela_electron_qf"]) - (item["sent_potential"] - item["sent_electron_qf"]))
        if item.get("sent_electron_density", 0.0) > 0.0 and item.get("vela_electron_density", 0.0) > 0.0:
            logn.append(math.log10(item["vela_electron_density"] / item["sent_electron_density"]))
    return {
        "median_delta_psi_minus_phin_V": median(dpsi_phin),
        "median_log10_n_vela_over_sentaurus": median(logn),
        "node_count": len(by_node),
    }


def point_summary(case_dir: Path, tag: str) -> dict[str, Any]:
    rows = read_rows(case_dir / f"coarse_qf_case_{tag}.csv")
    row = next((r for r in rows if close(fnum(r.get("bias_V")), -18.0)), rows[-1] if rows else {})
    return {
        "converged": row.get("converged", ""),
        "iterations": row.get("iterations", ""),
        "newton_iterations": row.get("newton_iterations", ""),
        "failure_reason": row.get("failure_reason", ""),
        "current_total": row.get("current_total", ""),
    }


def newton_summary(case_dir: Path, tag: str) -> dict[str, Any]:
    path = case_dir / f"newton_history_case_{tag}.csv"
    if not path.exists():
        return {"history_rows": 0, "history_max_iteration": ""}
    rows = [r for r in read_rows(path) if close(fnum(r.get("bias_V")), -18.0)]
    iters = [int(float(r["iteration"])) for r in rows if r.get("iteration")]
    return {
        "history_rows": len(rows),
        "history_max_iteration": max(iters) if iters else "",
    }


def edge_id(row: dict[str, str]) -> tuple[int, int, int]:
    return (int(row["edge_id"]), int(row["node0"]), int(row["node1"]))


def support_flux_metrics(row: dict[str, str]) -> tuple[float, float, float]:
    vela_flux = abs(fnum(row.get("electron_flux_abs")) or 0.0) + abs(fnum(row.get("hole_flux_abs")) or 0.0)
    q = 1.602176634e-19
    sent_flux = (
        abs(fnum(row.get("sent_e_current")) or 0.0) * 1.0e4 / q
        + abs(fnum(row.get("sent_h_current")) or 0.0) * 1.0e4 / q
    )
    ratio = math.inf if sent_flux == 0.0 and vela_flux > 0.0 else (vela_flux / sent_flux if sent_flux else math.nan)
    return vela_flux, sent_flux, ratio


def state_by_node(node_csv: Path) -> dict[int, dict[str, float]]:
    out: dict[int, dict[str, float]] = {}
    for row in read_rows(node_csv):
        if not close(fnum(row.get("bias_V")), -18.0):
            continue
        q = row.get("quantity", "").strip()
        if q not in {"potential", "electron_qf", "electron_density"}:
            continue
        node = int(row["node_id"])
        item = out.setdefault(node, {})
        vv = fnum(row.get("vela_value_scaled_to_sentaurus_units"))
        sv = fnum(row.get("sentaurus_value"))
        if vv is not None:
            item[f"vela_{q}"] = vv
        if sv is not None:
            item[f"sent_{q}"] = sv
    return out


def quantiles(values: list[float], probs: list[float]) -> dict[str, float | None]:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return {f"p{int(p * 100):02d}": None for p in probs}
    result: dict[str, float | None] = {}
    for p in probs:
        idx = min(len(clean) - 1, max(0, round(p * (len(clean) - 1))))
        result[f"p{int(p * 100):02d}"] = clean[idx]
    return result


def edge_step_analysis(case_dir: Path, tag: str, node_csv: Path, support_csv: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    states = state_by_node(node_csv)
    rows = [r for r in read_rows(support_csv) if close(fnum(r.get("bias_V")), -18.0)]
    all_edges: list[dict[str, Any]] = []
    for row in rows:
        n0 = int(row["node0"])
        n1 = int(row["node1"])
        if n0 not in states or n1 not in states:
            continue
        phin0 = states[n0].get("vela_electron_qf")
        phin1 = states[n1].get("vela_electron_qf")
        psi0 = states[n0].get("vela_potential")
        psi1 = states[n1].get("vela_potential")
        if phin0 is None or phin1 is None:
            continue
        vela_flux, sent_flux, flux_ratio = support_flux_metrics(row)
        e_alpha = fnum(row.get("electron_alpha_m_inv")) or 0.0
        h_alpha = fnum(row.get("hole_alpha_m_inv")) or 0.0
        e_qf = fnum(row.get("electron_qf_field_V_m"))
        h_qf = fnum(row.get("hole_qf_field_V_m"))
        electric = fnum(row.get("electric_field_V_m"))
        all_edges.append({
            "case": tag,
            "edge_id": row["edge_id"],
            "node0": n0,
            "node1": n1,
            "edge_class": row.get("edge_class", ""),
            "abs_dphin_V": abs(phin1 - phin0),
            "abs_dpsi_V": None if psi0 is None or psi1 is None else abs(psi1 - psi0),
            "vela_flux_m2_s": vela_flux,
            "sentaurus_flux_m2_s": sent_flux,
            "vela_over_sentaurus_flux": flux_ratio,
            "combined_alpha_m_inv": e_alpha + h_alpha,
            "electron_alpha_m_inv": e_alpha,
            "hole_alpha_m_inv": h_alpha,
            "electron_qf_field_V_m": e_qf,
            "hole_qf_field_V_m": h_qf,
            "electric_field_V_m": electric,
            "electron_qf_over_electric": (e_qf / electric if e_qf is not None and electric not in (None, 0.0) else None),
            "hole_qf_over_electric": (h_qf / electric if h_qf is not None and electric not in (None, 0.0) else None),
            "minimum_field_or_exp_zero_alpha": int((e_alpha + h_alpha) == 0.0),
        })
    dphins = [float(r["abs_dphin_V"]) for r in all_edges]
    qs = quantiles(dphins, [0.5, 0.9, 0.95, 0.99])
    for r in all_edges:
        r["is_flux_gt_1e3"] = int((r["vela_over_sentaurus_flux"] or 0.0) > 1.0e3)
        r["all_edge_abs_dphin_p50_V"] = qs["p50"]
        r["all_edge_abs_dphin_p90_V"] = qs["p90"]
        r["all_edge_abs_dphin_p95_V"] = qs["p95"]
        r["all_edge_abs_dphin_p99_V"] = qs["p99"]
    hot = sorted(
        [r for r in all_edges if r["is_flux_gt_1e3"]],
        key=lambda r: (-(r["vela_over_sentaurus_flux"] if math.isfinite(r["vela_over_sentaurus_flux"]) else 1.0e300), -r["abs_dphin_V"]),
    )
    return all_edges, hot


def plot_dphin_distribution(all_edges: list[dict[str, Any]], hot_edges: list[dict[str, Any]], out_path: Path) -> str | None:
    try:
        os.environ.setdefault("MPLCONFIGDIR", str(out_path.parent / ".matplotlib"))
        Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
        import matplotlib.pyplot as plt
    except Exception as exc:
        return str(exc)
    all_vals = [float(r["abs_dphin_V"]) for r in all_edges]
    hot_vals = [float(r["abs_dphin_V"]) for r in hot_edges]
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    bins = 30
    ax.hist(all_vals, bins=bins, color="#3A6EA5", alpha=0.55, label=f"all active edges (n={len(all_vals)})")
    if hot_vals:
        ax.hist(hot_vals, bins=min(bins, max(5, len(hot_vals))), color="#B55A30", alpha=0.75, label=f"flux ratio >1e3 (n={len(hot_vals)})")
    ax.set_xlabel("|phin_j - phin_i| on Vela edge (V)")
    ax.set_ylabel("edge count")
    ax.set_title("-18 V active-edge electron quasi-Fermi step distribution")
    ax.grid(axis="y", linewidth=0.45, alpha=0.35)
    ax.legend()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return None


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    baseline_all: list[dict[str, Any]] = []
    baseline_hot: list[dict[str, Any]] = []
    for tag, label, warm, cap, line_search, damping in CASES:
        case_dir = OUT_DIR / f"case_{tag}"
        if case_dir.exists():
            shutil.rmtree(case_dir)
        case_dir.mkdir(parents=True)
        config_path = configure_case(case_dir, tag, warm, cap, line_search, damping)
        returncode, _, _ = run_case(case_dir, config_path)
        node_csv = support_csv = None
        branch = {"median_delta_psi_minus_phin_V": None, "median_log10_n_vela_over_sentaurus": None, "node_count": 0}
        if (case_dir / f"sg_avalanche_edges_case_{tag}.csv").exists():
            node_csv, support_csv = build_compare_outputs(case_dir, tag)
            branch = branch_summary_from_node_csv(node_csv)
        point = point_summary(case_dir, tag) if (case_dir / f"coarse_qf_case_{tag}.csv").exists() else {}
        hist = newton_summary(case_dir, tag)
        if tag == "A" and node_csv is not None and support_csv is not None:
            baseline_all, baseline_hot = edge_step_analysis(case_dir, tag, node_csv, support_csv)
        summary_rows.append({
            "case": tag,
            "label": label,
            "runner_returncode": returncode,
            "warm_start": warm,
            "quasi_fermi_update_limit_V": cap,
            "line_search": line_search,
            "damping_factor": damping,
            **point,
            **hist,
            **branch,
        })

    write_rows(OUT_DIR / "qf_cap_warmstart_summary.csv", summary_rows, [
        "case", "label", "runner_returncode", "warm_start", "quasi_fermi_update_limit_V",
        "line_search", "damping_factor", "converged", "iterations", "newton_iterations",
        "history_rows", "history_max_iteration", "failure_reason", "current_total",
        "median_delta_psi_minus_phin_V", "median_log10_n_vela_over_sentaurus", "node_count",
    ])
    if baseline_all:
        fields = [
            "case", "edge_id", "node0", "node1", "edge_class", "abs_dphin_V", "abs_dpsi_V",
            "vela_flux_m2_s", "sentaurus_flux_m2_s", "vela_over_sentaurus_flux",
            "combined_alpha_m_inv", "electron_alpha_m_inv", "hole_alpha_m_inv",
            "electron_qf_field_V_m", "hole_qf_field_V_m", "electric_field_V_m",
            "electron_qf_over_electric", "hole_qf_over_electric",
            "minimum_field_or_exp_zero_alpha", "is_flux_gt_1e3",
            "all_edge_abs_dphin_p50_V", "all_edge_abs_dphin_p90_V",
            "all_edge_abs_dphin_p95_V", "all_edge_abs_dphin_p99_V",
        ]
        write_rows(OUT_DIR / "baseline_active_edge_dphin_all.csv", baseline_all, fields)
        write_rows(OUT_DIR / "baseline_flux_gt_1e3_edge_dphin.csv", baseline_hot, fields)
        plot_error = plot_dphin_distribution(
            baseline_all,
            baseline_hot,
            OUT_DIR / "baseline_active_edge_dphin_distribution.png",
        )
        if plot_error:
            (OUT_DIR / "plot_error.txt").write_text(plot_error + "\n", encoding="utf-8")
    print(json.dumps({"out_dir": str(OUT_DIR), "cases": len(summary_rows), "baseline_hot_edges": len(baseline_hot)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())




