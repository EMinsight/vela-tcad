#!/usr/bin/env python3
"""Replay exported Sentaurus TransportModels states through Vela formulas.

The replay is immutable: Sentaurus electrostatic/quasi-Fermi/quantum fields are
loaded as an external state and Vela evaluates carrier density, mobility, SG
edge flux, and SRH continuity terms without taking a nonlinear solve step.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
TRANSITION_MANIFEST = REF / "sentaurus_vm_runs/idvg_spatial_oracle_20260821/spatial_oracle_manifest.json"
REMAINING_MANIFEST = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/remaining_spatial_oracles_manifest.json"
BASELINE = REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23/runs/dg"
IDVG_CONFIG = BASELINE / "03_dg_idvg_curve.json"
IDVD_CONFIG = BASELINE / "05_dg_idvd_curve.json"
DEFAULT_RUNNER = REPO / "build/vela_example_runner.exe"
OUTPUT = REF / "reports/transportmodels_sentaurus_formula_replay_20260823"
REPORT_JSON = REPO / "docs/validation/transportmodels_sentaurus_formula_replay_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_sentaurus_formula_replay_2026-08-23.md"
Q = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def error_stats(values: Iterable[float]) -> dict[str, float]:
    data = [float(value) for value in values if math.isfinite(float(value))]
    return {
        "count": len(data),
        "median": percentile(data, 0.5),
        "p95": percentile(data, 0.95),
        "maximum": max(data, default=math.nan),
    }


def field(export_dir: Path, name: str, region: int, components: int = 1) -> dict[int, Any]:
    path = export_dir / "fields" / f"{name}_region{region}.csv"
    rows = read_csv(path)
    if components == 1:
        return {int(row["node_id"]): float(row["component0"]) for row in rows}
    return {
        int(row["node_id"]): tuple(float(row[f"component{i}"]) for i in range(components))
        for row in rows
    }


def slug(kind: str, value: float) -> str:
    return f"{kind}_{value:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")


def cases() -> list[dict[str, Any]]:
    transition = json.loads(TRANSITION_MANIFEST.read_text(encoding="utf-8"))
    remaining = json.loads(REMAINING_MANIFEST.read_text(encoding="utf-8"))
    result: list[dict[str, Any]] = []
    for state in transition["states"]:
        bias = float(state["gate_bias_V"])
        result.append(
            {
                "group": "dg_idvg_transition",
                "bias_kind": "gate",
                "bias_V": bias,
                "gate_bias_V": bias,
                "drain_bias_V": 1.1,
                "export_dir": Path(state["export_dir"]),
            }
        )
    for state in remaining["idvd_states"]:
        bias = float(state["drain_bias_V"])
        result.append(
            {
                "group": "dg_idvd",
                "bias_kind": "drain",
                "bias_V": bias,
                "gate_bias_V": 1.0,
                "drain_bias_V": bias,
                "export_dir": Path(state["export_dir"]),
            }
        )
    for state in remaining["dg_deep_off_states"]:
        bias = float(state["gate_bias_V"])
        result.append(
            {
                "group": "dg_idvg_deep_off",
                "bias_kind": "gate",
                "bias_V": bias,
                "gate_bias_V": bias,
                "drain_bias_V": 1.1,
                "export_dir": Path(state["export_dir"]),
            }
        )
    return result


def runner_environment() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    return env


def convert_state(export_dir: Path, output: Path) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/sentaurus_fields_to_restart.py"),
            "--export-dir",
            str(export_dir),
            "--output",
            str(output),
            "--preserve-insulator-quantum-potential",
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)


def feedback_fields(state_file: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv(state_file)
    mapping = {
        "eQuasiFermiPotential": "phin",
        "hQuasiFermiPotential": "phip",
        "eDensity_m3": "electrons_m3",
        "hDensity_m3": "holes_m3",
    }
    for name, column in mapping.items():
        density = name.endswith("Density_m3")
        write_csv(
            output_dir / f"{name}_region0.csv",
            [
                {
                    "node_id": int(row["node_id"]),
                    "component0": max(float(row[column]), 1.0) if density else float(row[column]),
                }
                for row in rows
            ],
        )
    return output_dir


def make_probe_config(
    case: dict[str, Any], simulation_type: str, state_file: Path, output_csv: Path,
    feedback_dir: Path | None = None,
) -> dict[str, Any]:
    source = IDVD_CONFIG if case["group"] == "dg_idvd" else IDVG_CONFIG
    config = json.loads(source.read_text(encoding="utf-8"))
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state_file.resolve())
    config["output_csv"] = str(output_csv.resolve())
    config.pop("sweep", None)
    config.pop("log_file", None)
    for contact in config["contacts"]:
        name = contact["name"].lower()
        if name == "gate":
            contact["bias"] = case["gate_bias_V"]
        elif name == "drain":
            contact["bias"] = case["drain_bias_V"]
    if feedback_dir is not None:
        config["feedback_state_fields_dir"] = str(feedback_dir.resolve())
    config["_comment"] = (
        "Immutable Sentaurus state replay through Vela production formulas; "
        f"group={case['group']}, {case['bias_kind']}={case['bias_V']} V"
    )
    return config


def run_probe(
    runner: Path, case: dict[str, Any], run_dir: Path, simulation_type: str,
    state_file: Path, feedback_dir: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    label = simulation_type.removesuffix("_probe")
    output_csv = run_dir / f"{label}.csv"
    config_path = run_dir / f"{label}.json"
    config_path.write_text(
        json.dumps(
            make_probe_config(case, simulation_type, state_file, output_csv, feedback_dir),
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [str(runner), "--config", str(config_path), "--log", str(run_dir / f"{label}.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (run_dir / f"{label}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_dir / f"{label}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"{simulation_type} failed: {completed.stderr or completed.stdout}")
    status = json.loads(completed.stdout.strip().splitlines()[-1])
    return output_csv, status


def endpoint_density_metrics(export_dir: Path, edge_rows: list[dict[str, str]]) -> dict[str, Any]:
    sent_n = field(export_dir, "eDensity", 3)
    sent_p = field(export_dir, "hDensity", 3)
    vela_n: dict[int, list[float]] = {}
    vela_p: dict[int, list[float]] = {}
    for row in edge_rows:
        for endpoint in (0, 1):
            node = int(row[f"node{endpoint}"])
            if node in sent_n:
                vela_n.setdefault(node, []).append(float(row[f"electron_density{endpoint}_m3"]) / 1.0e6)
                vela_p.setdefault(node, []).append(float(row[f"hole_density{endpoint}_m3"]) / 1.0e6)
    n_errors = [
        abs(math.log10(max(statistics.median(vela_n[node]), 1.0)) - math.log10(max(sent_n[node], 1.0)))
        for node in vela_n
    ]
    p_errors = [
        abs(math.log10(max(statistics.median(vela_p[node]), 1.0)) - math.log10(max(sent_p[node], 1.0)))
        for node in vela_p
    ]
    return {"electron_density_abs_error_dex": error_stats(n_errors), "hole_density_abs_error_dex": error_stats(p_errors)}


def edge_transport_metrics(
    export_dir: Path,
    sg_rows: list[dict[str, str]],
    mobility_rows: list[dict[str, str]],
    comparison_csv: Path,
) -> dict[str, Any]:
    sent_mobility = field(export_dir, "eMobility", 3)
    sent_drive = field(export_dir, "eEparallel", 3)
    sent_current = field(export_dir, "eCurrentDensity", 3, 2)
    mobility_errors: list[float] = []
    mobility_error_by_edge: dict[int, float] = {}
    drive_pairs: list[tuple[float, float]] = []
    drive_error_by_edge: dict[int, float] = {}
    mobility_by_edge: dict[int, dict[str, float]] = {}
    for row in mobility_rows:
        edge_id = int(row["edge_id"])
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 not in sent_mobility or n1 not in sent_mobility:
            continue
        sent_mu = 0.5 * (sent_mobility[n0] + sent_mobility[n1])
        vela_mu = float(row["electron_final_mobility_m2_V_s"]) * 1.0e4
        mobility_error = abs(math.log10(max(vela_mu, 1.0e-30)) - math.log10(max(sent_mu, 1.0e-30)))
        mobility_errors.append(mobility_error)
        mobility_error_by_edge[edge_id] = mobility_error
        sent_e = 0.5 * (abs(sent_drive[n0]) + abs(sent_drive[n1]))
        vela_e = float(row["electron_mobility_field_V_m"]) / 100.0
        drive_pairs.append((sent_e, vela_e))
        if sent_e > 0.0 and vela_e > 0.0:
            drive_error_by_edge[edge_id] = abs(math.log10(vela_e) - math.log10(sent_e))
        mobility_by_edge[edge_id] = {
            "sentaurus_eMobility_cm2_V_s": sent_mu,
            "vela_eMobility_cm2_V_s": vela_mu,
            "sentaurus_eEparallel_V_cm": sent_e,
            "vela_eMobilityDrive_V_cm": vela_e,
        }

    drive_floor = max((max(pair) for pair in drive_pairs), default=0.0) * 1.0e-10
    drive_errors = [
        abs(math.log10(vela) - math.log10(sent))
        for sent, vela in drive_pairs
        if sent > drive_floor and vela > drive_floor
    ]

    current_pairs: list[tuple[float, float]] = []
    edge_rows: list[dict[str, Any]] = []
    for row in sg_rows:
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 not in sent_current or n1 not in sent_current:
            continue
        dx = float(row["x1"]) - float(row["x0"])
        dy = float(row["y1"]) - float(row["y0"])
        length = math.hypot(dx, dy)
        if length == 0.0:
            continue
        tx, ty = dx / length, dy / length
        jx = 0.5 * (sent_current[n0][0] + sent_current[n1][0])
        jy = 0.5 * (sent_current[n0][1] + sent_current[n1][1])
        sent_line_A_m = abs((jx * tx + jy * ty) * 1.0e4 * float(row["couple_m"]))
        vela_line_A_m = abs(Q * float(row["electron_particle_line_flux_per_m_s"]))
        current_pairs.append((sent_line_A_m, vela_line_A_m))
        edge_id = int(row["edge_id"])
        extra = mobility_by_edge.get(edge_id, {})
        edge_rows.append(
            {
                "edge_id": edge_id,
                "node0": n0,
                "node1": n1,
                "x0_um": float(row["x0"]) * 1.0e6,
                "y0_um": float(row["y0"]) * 1.0e6,
                "x1_um": float(row["x1"]) * 1.0e6,
                "y1_um": float(row["y1"]) * 1.0e6,
                "sentaurus_eLineCurrent_A_m": sent_line_A_m,
                "vela_eSgLineCurrent_A_m": vela_line_A_m,
                "line_current_abs_error_dex": (
                    abs(math.log10(vela_line_A_m) - math.log10(sent_line_A_m))
                    if sent_line_A_m > 0.0 and vela_line_A_m > 0.0 else math.nan
                ),
                **extra,
                "mobility_abs_error_dex": mobility_error_by_edge.get(edge_id, math.nan),
                "drive_abs_error_dex": drive_error_by_edge.get(edge_id, math.nan),
            }
        )
    write_csv(comparison_csv, edge_rows)
    current_floor = max((max(pair) for pair in current_pairs), default=0.0) * 1.0e-10
    current_errors = [
        abs(math.log10(vela) - math.log10(sent))
        for sent, vela in current_pairs
        if sent > current_floor and vela > current_floor
    ]
    max_sent_current = max((pair[0] for pair in current_pairs), default=0.0)
    active_edge_ids = {
        int(row["edge_id"])
        for row in edge_rows
        if float(row["sentaurus_eLineCurrent_A_m"]) > max_sent_current * 1.0e-3
    }
    active_current_errors = [
        float(row["line_current_abs_error_dex"])
        for row in edge_rows
        if int(row["edge_id"]) in active_edge_ids
        and math.isfinite(float(row["line_current_abs_error_dex"]))
    ]
    active_mobility_errors = [
        mobility_error_by_edge[edge] for edge in active_edge_ids if edge in mobility_error_by_edge
    ]
    active_drive_errors = [
        drive_error_by_edge[edge] for edge in active_edge_ids if edge in drive_error_by_edge
    ]
    return {
        "electron_mobility_abs_error_dex": error_stats(mobility_errors),
        "electron_drive_abs_error_dex": error_stats(drive_errors),
        "electron_sg_line_current_abs_error_dex": error_stats(current_errors),
        "current_carrying_edges": {
            "definition": "Sentaurus projected electron line current exceeds 1e-3 of its maximum",
            "edge_count": len(active_edge_ids),
            "electron_mobility_abs_error_dex": error_stats(active_mobility_errors),
            "electron_drive_abs_error_dex": error_stats(active_drive_errors),
            "electron_sg_line_current_abs_error_dex": error_stats(active_current_errors),
        },
        "sg_current_comparison_note": "Sentaurus nodal J projected on the primal edge and integrated over the Vela dual couple; use as a spatial diagnostic, not an exact discretization identity.",
    }


def srh_shape_metrics(export_dir: Path, term_rows: list[dict[str, str]]) -> dict[str, Any]:
    sent_srh = field(export_dir, "srhRecombination", 3)
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(export_dir / "nodes.csv")
    }
    area = {node: 0.0 for node in sent_srh}
    for element in read_csv(export_dir / "elements.csv"):
        if element["region"] != "R.Substrate":
            continue
        ids = tuple(int(element[f"node{i}"]) for i in range(3))
        (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in ids)
        triangle_area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        for node in ids:
            area[node] += triangle_area / 3.0
    sent_weight = {node: abs(sent_srh[node] * area[node]) for node in sent_srh}
    vela_weight = {
        int(row["node_id"]): abs(float(row["electron_recombination"]))
        for row in term_rows
        if int(row["node_id"]) in sent_srh
    }
    sent_total = sum(sent_weight.values())
    vela_total = sum(vela_weight.values())
    nodes_common = set(sent_weight) & set(vela_weight)
    tv = 0.5 * sum(
        abs(sent_weight[node] / max(sent_total, 1.0e-300) - vela_weight[node] / max(vela_total, 1.0e-300))
        for node in nodes_common
    )
    sent_top = {node for node, _ in sorted(sent_weight.items(), key=lambda item: item[1], reverse=True)[:50]}
    vela_top = {node for node, _ in sorted(vela_weight.items(), key=lambda item: item[1], reverse=True)[:50]}
    return {
        "normalized_total_variation_distance": tv,
        "top50_node_overlap": len(sent_top & vela_top) / 50.0,
        "sentaurus_signed_area_weighted_sum_cm-1_s-1": sum(sent_srh[node] * area[node] for node in sent_srh),
        "vela_signed_scaled_integral": sum(
            float(row["electron_recombination"])
            for row in term_rows
            if int(row["node_id"]) in sent_srh
        ),
    }


def feedback_metrics(rows: list[dict[str, str]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["variant"], []).append(row)
    return {
        variant: {
            "electron_recombination_sum": sum(float(row["electron_recombination"]) for row in values),
            "electron_flux_sum": sum(float(row["electron_flux"]) for row in values),
            "electron_residual_l2": math.sqrt(sum(float(row["electron_residual"]) ** 2 for row in values)),
        }
        for variant, values in grouped.items()
    }


def execute_case(runner: Path, case: dict[str, Any], replay: bool) -> dict[str, Any]:
    run_dir = OUTPUT / case["group"] / slug(case["bias_kind"], case["bias_V"])
    run_dir.mkdir(parents=True, exist_ok=True)
    state_file = run_dir / "sentaurus_state_for_vela.csv"
    previous_path = run_dir / "formula_replay_summary.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8")) if previous_path.is_file() else {}
    if replay:
        convert_state(case["export_dir"], state_file)
        sg_csv, sg_status = run_probe(runner, case, run_dir, "sg_edge_flux_probe", state_file)
        mobility_csv, mobility_status = run_probe(runner, case, run_dir, "edge_mobility_probe", state_file)
        term_csv, term_status = run_probe(runner, case, run_dir, "newton_carrier_term_probe", state_file)
    else:
        sg_csv = run_dir / "sg_edge_flux.csv"
        mobility_csv = run_dir / "edge_mobility.csv"
        term_csv = run_dir / "newton_carrier_term.csv"
        sg_status = previous.get("statuses", {}).get("sg", {})
        mobility_status = previous.get("statuses", {}).get("mobility", {})
        term_status = previous.get("statuses", {}).get("carrier_terms", {})
        for path in (state_file, sg_csv, mobility_csv, term_csv):
            if not path.is_file():
                raise FileNotFoundError(f"Analyze-only input is missing: {path}")
    sg_rows, mobility_rows, term_rows = read_csv(sg_csv), read_csv(mobility_csv), read_csv(term_csv)
    result = {
        "group": case["group"],
        "bias_kind": case["bias_kind"],
        "bias_V": case["bias_V"],
        "gate_bias_V": case["gate_bias_V"],
        "drain_bias_V": case["drain_bias_V"],
        "density": endpoint_density_metrics(case["export_dir"], sg_rows),
        "transport": edge_transport_metrics(
            case["export_dir"], sg_rows, mobility_rows,
            run_dir / "sentaurus_vela_edge_formula_comparison.csv",
        ),
        "srh": srh_shape_metrics(case["export_dir"], term_rows),
        "statuses": {"sg": sg_status, "mobility": mobility_status, "carrier_terms": term_status},
        "artifacts": {
            "sentaurus_state": str(state_file.resolve()),
            "sg_edges": str(sg_csv.resolve()),
            "edge_mobility": str(mobility_csv.resolve()),
            "edge_formula_comparison": str(
                (run_dir / "sentaurus_vela_edge_formula_comparison.csv").resolve()
            ),
            "carrier_terms": str(term_csv.resolve()),
        },
    }
    if case["group"] == "dg_idvg_deep_off":
        fields = feedback_fields(state_file, run_dir / "feedback_fields")
        if replay:
            feedback_csv, feedback_status = run_probe(
                runner, case, run_dir, "newton_feedback_substitution_probe", state_file, fields
            )
        else:
            feedback_csv = run_dir / "newton_feedback_substitution.csv"
            feedback_status = previous.get("statuses", {}).get("feedback", {})
        result["feedback_substitution"] = feedback_metrics(read_csv(feedback_csv))
        result["statuses"]["feedback"] = feedback_status
        result["artifacts"]["feedback"] = str(feedback_csv.resolve())
    summary_path = run_dir / "formula_replay_summary.json"
    summary_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    result["artifacts"]["summary"] = str(summary_path.resolve())
    return result


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels Sentaurus fixed-state formula replay",
        "",
        "Sentaurus 2022 fields are inserted into Vela production density, mobility, SG flux, and SRH operators without a nonlinear update.",
        "",
        "| Group | Bias | n p95 (dex) | mobility p95 (all/active) | SG J p95 (all/active) | SRH shape TV |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in report["cases"]:
        lines.append(
            f"| {row['group']} | {row['bias_V']:.2f} V | "
            f"{row['density']['electron_density_abs_error_dex']['p95']:.4g} | "
            f"{row['transport']['electron_mobility_abs_error_dex']['p95']:.4g} / "
            f"{row['transport']['current_carrying_edges']['electron_mobility_abs_error_dex']['p95']:.4g} | "
            f"{row['transport']['electron_sg_line_current_abs_error_dex']['p95']:.4g} / "
            f"{row['transport']['current_carrying_edges']['electron_sg_line_current_abs_error_dex']['p95']:.4g} | "
            f"{row['srh']['normalized_total_variation_distance']:.4g} |"
        )
    lines.extend(
        [
            "",
            "The SG comparison projects Sentaurus nodal current density onto Vela primal edges and integrates over the Vela dual couple. It is a localization diagnostic rather than an assertion that the two discretizations are identical.",
            "",
            f"Raw artifact directory: `{report['output_directory']}`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for row in report["cases"]:
            for key in ("sentaurus_state", "sg_edges", "edge_mobility", "carrier_terms", "summary"):
                if not Path(row["artifacts"][key]).is_file():
                    raise RuntimeError(f"Missing artifact {key}: {row['artifacts'][key]}")
        print("TransportModels Sentaurus formula replay check: PASS")
        return 0

    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = []
    for index, case in enumerate(cases(), start=1):
        print(f"[{index}/12] {case['group']} {case['bias_kind']}={case['bias_V']:.2f} V", flush=True)
        results.append(execute_case(args.runner.resolve(), case, not args.analyze_only))
    report = {
        "schema": "vela.transportmodels.sentaurus_formula_replay.v1",
        "as_of": "2026-08-23",
        "status": "complete",
        "runner": str(args.runner.resolve()),
        "runner_sha256": sha256(args.runner.resolve()),
        "case_count": len(results),
        "cases": results,
        "output_directory": str(OUTPUT.resolve()),
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": "complete", "cases": len(results), "report": str(REPORT_JSON)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
