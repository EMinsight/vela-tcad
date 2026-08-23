#!/usr/bin/env python3
"""Audit the deep-off Poisson line search and local carrier-row coupling."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any


ALPHAS = [2.0 ** (-i) for i in range(13)]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def replace_path_prefix(value: Any, old: str, new: str) -> Any:
    if isinstance(value, dict):
        return {key: replace_path_prefix(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_path_prefix(item, old, new) for item in value]
    if isinstance(value, str):
        return value.replace(old, new)
    return value


def run_probe(runner: Path, config: dict[str, Any], path: Path) -> dict[str, Any]:
    write_json(path, config)
    completed = subprocess.run(
        [str(runner), "--config", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    status: dict[str, Any] = {}
    for line in reversed(completed.stdout.splitlines()):
        line = line.strip()
        if line.startswith("{"):
            status = json.loads(line)
            break
    status["returncode"] = completed.returncode
    status["stdout"] = completed.stdout
    status["stderr"] = completed.stderr
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed ({completed.returncode}) for {path}:\n{completed.stderr}"
        )
    return status


def probe_config(
    base: dict[str, Any],
    simulation_type: str,
    state_file: Path,
    output_csv: Path,
    gate_bias: float = -0.68,
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state_file.resolve())
    config["output_csv"] = str(output_csv.resolve())
    config.pop("sweep", None)
    for contact in config.get("contacts", []):
        if str(contact.get("name", "")).lower() == "gate":
            contact["bias"] = gate_bias
    return config


def write_candidate_state(
    base_rows: list[dict[str, str]],
    step_rows: list[dict[str, str]],
    alpha: float,
    path: Path,
) -> None:
    step_by_node = {int(row["node_id"]): row for row in step_rows}
    fieldnames = list(base_rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for source in base_rows:
            row = dict(source)
            node = int(row["node_id"])
            step = step_by_node[node]
            row["psi"] = format(
                float(row["psi"]) + alpha * float(step["delta_psi_V"]), ".17g"
            )
            electron_increment = (
                float(row["electron_qf_increment_V"])
                + alpha * float(step["delta_phin_V"])
            )
            hole_increment = (
                float(row["hole_qf_increment_V"])
                + alpha * float(step["delta_phip_V"])
            )
            electron_reference = float(row["electron_qf_reference_V"])
            hole_reference = float(row["hole_qf_reference_V"])
            row["electron_qf_increment_V"] = format(electron_increment, ".17g")
            row["hole_qf_increment_V"] = format(hole_increment, ".17g")
            row["phin"] = format(electron_reference + electron_increment, ".17g")
            row["phip"] = format(hole_reference + hole_increment, ".17g")
            writer.writerow(row)


def l2(values: list[float]) -> float:
    return math.sqrt(math.fsum(value * value for value in values))


def block_norms(rows: list[dict[str, str]]) -> dict[str, float]:
    psi = l2([float(row["psi_residual"]) for row in rows])
    phin = l2([float(row["phin_residual"]) for row in rows])
    phip = l2([float(row["phip_residual"]) for row in rows])
    return {
        "psi": psi,
        "phin": phin,
        "phip": phip,
        "combined": math.sqrt(psi * psi + phin * phin + phip * phip),
    }


def carrier_metrics(
    rows: list[dict[str, str]],
    eps_row: float,
    scale_floor: float,
    min_source_scale: float,
    min_source_fraction: float,
) -> tuple[list[dict[str, Any]], dict[int, float]]:
    violations: list[dict[str, Any]] = []
    node_max_ratio: dict[int, float] = {}
    for row in rows:
        node = int(row["node_id"])
        for carrier in ("electron", "hole"):
            residual = float(row[f"{carrier}_residual"])
            flux = float(row[f"{carrier}_flux"])
            flux_abs_sum = float(row[f"{carrier}_flux_abs_sum"])
            recombination = float(row[f"{carrier}_recombination"])
            impact = float(row[f"{carrier}_impact"])
            source_scale = max(abs(recombination), abs(impact))
            scale = max(abs(flux_abs_sum), source_scale, scale_floor)
            source_qualified = (
                scale > 0.0
                and source_scale >= min_source_scale
                and source_scale >= min_source_fraction * scale
            )
            ratio = abs(residual) / scale if scale > 0.0 else 0.0
            if source_qualified:
                node_max_ratio[node] = max(node_max_ratio.get(node, 0.0), ratio)
            if source_qualified and ratio > eps_row:
                violations.append(
                    {
                        "node_id": node,
                        "carrier": carrier,
                        "ratio": ratio,
                        "residual": residual,
                        "scale": scale,
                        "flux": flux,
                        "flux_abs_sum": flux_abs_sum,
                        "recombination": recombination,
                        "impact": impact,
                        "x": float(row["x"]),
                        "y": float(row["y"]),
                    }
                )
    violations.sort(key=lambda row: row["ratio"], reverse=True)
    return violations, node_max_ratio


def pearson_log_correlation(
    residual_rows: list[dict[str, str]], node_ratios: dict[int, float]
) -> float | None:
    pairs: list[tuple[float, float]] = []
    for row in residual_rows:
        node = int(row["node_id"])
        ratio = node_ratios.get(node)
        poisson = abs(float(row["psi_residual"]))
        if ratio is not None and ratio > 0.0 and poisson > 0.0:
            pairs.append((math.log10(poisson), math.log10(ratio)))
    if len(pairs) < 2:
        return None
    mean_x = math.fsum(x for x, _ in pairs) / len(pairs)
    mean_y = math.fsum(y for _, y in pairs) / len(pairs)
    covariance = math.fsum((x - mean_x) * (y - mean_y) for x, y in pairs)
    variance_x = math.fsum((x - mean_x) ** 2 for x, _ in pairs)
    variance_y = math.fsum((y - mean_y) ** 2 for _, y in pairs)
    denominator = math.sqrt(variance_x * variance_y)
    return covariance / denominator if denominator > 0.0 else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--runner", type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve()
    runner = (args.runner or repo / "build-qf-ab" / "vela_example_runner.exe").resolve()
    source = (
        repo
        / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
        / "idvg_deep_off_precision_20260822/contact_basin_full_noglobal_20260822"
    )
    output = source.parent / "poisson_carrier_coupling_audit_20260823"
    output.mkdir(parents=True, exist_ok=True)
    base_config = json.loads((source / "config.json").read_text(encoding="utf-8"))

    # Recreate the failure with the current executable.  The historical state
    # predates some reference-coordinate diagnostics, so a fresh in-memory
    # solve is the only reliable source for an exact line-search replay.
    formal_replay = output / "formal_replay_enforced_current"
    formal_replay.mkdir(parents=True, exist_ok=True)
    formal_config = replace_path_prefix(
        copy.deepcopy(base_config), str(source.resolve()), str(formal_replay.resolve())
    )
    formal_config["solver"]["carrier_row_convergence"]["mode"] = "enforce"
    formal_config_path = formal_replay / "config.json"
    write_json(formal_config_path, formal_config)
    formal_completed = subprocess.run(
        [str(runner), "--config", str(formal_config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    (formal_replay / "stdout.txt").write_text(
        formal_completed.stdout, encoding="utf-8"
    )
    (formal_replay / "stderr.txt").write_text(
        formal_completed.stderr, encoding="utf-8"
    )
    if formal_completed.returncode not in (0, 1):
        raise RuntimeError(
            f"formal replay failed ({formal_completed.returncode}):\n"
            f"{formal_completed.stderr}"
        )
    state_file = (
        formal_replay
        / "rejected_states/attempt_1_bias_m0p680000_final.csv"
    )
    if not state_file.exists():
        raise RuntimeError(f"formal replay did not write rejected state: {state_file}")
    formal_failures = json.loads(
        (formal_replay / "curve_newton_failure_diagnostics.json").read_text(
            encoding="utf-8"
        )
    )
    formal_failure = formal_failures[-1]
    formal_violation_rows = read_rows(formal_replay / "carrier_row_violations.csv")
    formal_event = formal_violation_rows[-1]["event"]
    formal_iteration = formal_violation_rows[-1]["iteration"]
    formal_violation_rows = [
        row
        for row in formal_violation_rows
        if row["event"] == formal_event and row["iteration"] == formal_iteration
    ]
    formal_history_rows = read_rows(formal_replay / "newton_history.csv")
    expected_formal_violations = next(
        int(row["carrier_row_violations"])
        for row in reversed(formal_history_rows)
        if row.get("carrier_row_violations", "").isdigit()
    )
    formal_violation_rows = formal_violation_rows[-expected_formal_violations:]
    base_config = formal_config

    full_step_csv = output / "full_step.csv"
    full_config = probe_config(base_config, "newton_step_probe", state_file, full_step_csv)
    full_status = run_probe(runner, full_config, output / "full_step_config.json")

    block_csv = output / "block_steps.csv"
    block_config = probe_config(base_config, "newton_block_step_probe", state_file, block_csv)
    block_config["block_modes"] = ["poisson_only", "carrier_only"]
    block_status = run_probe(runner, block_config, output / "block_steps_config.json")

    base_rows = read_rows(state_file)
    step_rows = read_rows(full_step_csv)
    carrier_cfg = base_config["solver"]["carrier_row_convergence"]
    eps_row = float(carrier_cfg.get("eps_row", 1.0e-5))
    scale_floor = float(carrier_cfg.get("scale_floor", 0.0))
    min_source_scale = float(carrier_cfg.get("min_source_scale", 0.0))
    min_source_fraction = float(carrier_cfg.get("min_source_scale_fraction", 0.0))

    candidate_rows: list[dict[str, Any]] = []
    baseline_residual_rows: list[dict[str, str]] = []
    baseline_ratios: dict[int, float] = {}
    baseline_violations: list[dict[str, Any]] = []
    for index, alpha in enumerate([0.0] + ALPHAS):
        label = "0" if alpha == 0.0 else f"{alpha:.12g}".replace(".", "p")
        candidate = output / f"candidate_alpha_{label}.csv"
        if alpha == 0.0:
            candidate = state_file
        else:
            write_candidate_state(base_rows, step_rows, alpha, candidate)

        residual_csv = output / f"residual_alpha_{label}.csv"
        residual_config = probe_config(
            base_config, "newton_residual_probe", candidate, residual_csv
        )
        run_probe(runner, residual_config, output / f"residual_alpha_{label}_config.json")
        residual_rows = read_rows(residual_csv)
        norms = block_norms(residual_rows)

        terms_csv = output / f"carrier_terms_alpha_{label}.csv"
        terms_config = probe_config(
            base_config, "newton_carrier_term_probe", candidate, terms_csv
        )
        terms_config["carrier_term_probe"] = {"solved_equation_terms": True}
        run_probe(runner, terms_config, output / f"carrier_terms_alpha_{label}_config.json")
        violations, node_ratios = carrier_metrics(
            read_rows(terms_csv), eps_row, scale_floor,
            min_source_scale, min_source_fraction
        )
        if alpha == 0.0:
            baseline_residual_rows = residual_rows
            baseline_ratios = node_ratios
            baseline_violations = violations
        node_1257_residual = next(
            row for row in residual_rows if int(row["node_id"]) == 1257
        )
        candidate_rows.append(
            {
                "alpha": alpha,
                **norms,
                "node_1257_psi_V": float(node_1257_residual["psi"]),
                "node_1257_psi_residual": float(
                    node_1257_residual["psi_residual"]
                ),
                "carrier_violations": len(violations),
                "carrier_max_ratio": violations[0]["ratio"] if violations else 0.0,
                "carrier_max_ratio_node": violations[0]["node_id"] if violations else -1,
                "carrier_max_ratio_type": violations[0]["carrier"] if violations else "",
            }
        )

    baseline_norm = candidate_rows[0]["combined"]
    for row in candidate_rows:
        row["relative_combined"] = row["combined"] / baseline_norm
        row["strict_decrease"] = row["combined"] < baseline_norm

    with (output / "candidate_merit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(candidate_rows[0].keys()))
        writer.writeheader()
        writer.writerows(candidate_rows)

    top_poisson = sorted(
        (
            {
                "node_id": int(row["node_id"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "poisson_residual": float(row["psi_residual"]),
            }
            for row in baseline_residual_rows
        ),
        key=lambda row: abs(row["poisson_residual"]),
        reverse=True,
    )
    violation_nodes = {row["node_id"] for row in baseline_violations}
    top10_nodes = {row["node_id"] for row in top_poisson[:10]}
    top50_nodes = {row["node_id"] for row in top_poisson[:50]}
    step_by_node = {int(row["node_id"]): row for row in step_rows}
    residual_by_node = {
        int(row["node_id"]): row for row in baseline_residual_rows
    }
    top_poisson_with_steps = []
    for row in top_poisson[:50]:
        step = step_by_node[row["node_id"]]
        base_psi = float(residual_by_node[row["node_id"]]["psi"])
        psi_ulp = math.ulp(base_psi)
        delta_psi = float(step["delta_psi_V"])
        top_poisson_with_steps.append(
            {
                **row,
                "psi_V": base_psi,
                "psi_ulp_V": psi_ulp,
                "delta_psi_V": delta_psi,
                "delta_psi_over_ulp": delta_psi / psi_ulp,
                "trial_poisson_residual": float(step["trial_psi_residual"]),
                "carrier_violation": row["node_id"] in violation_nodes,
                "carrier_max_ratio": baseline_ratios.get(row["node_id"], 0.0),
            }
        )
    with (output / "top_poisson_nodes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(top_poisson_with_steps[0].keys()))
        writer.writeheader()
        writer.writerows(top_poisson_with_steps)
    with (output / "carrier_violations_alpha_0.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(baseline_violations[0].keys()))
        writer.writeheader()
        writer.writerows(baseline_violations)

    summary = {
        "source_state": str(state_file.resolve()),
        "formal_replay_returncode": formal_completed.returncode,
        "formal_in_memory": {
            "failure_reason": formal_failure["failure_reason"],
            "block_residuals": formal_failure["block_residuals"],
            "carrier_violations": len(formal_violation_rows),
            "carrier_max_ratio": max(
                float(row["ratio"]) for row in formal_violation_rows
            ),
            "line_search_history": formal_failure["line_search_history"],
            "top10_poisson_carrier_overlap": len(
                {
                    int(row["node_id"])
                    for row in formal_failure["top_poisson_residual_nodes"][:10]
                }
                & {int(row["node_id"]) for row in formal_violation_rows}
            ),
        },
        "full_step": full_status,
        "block_steps": block_status,
        "candidate_merit": candidate_rows,
        "baseline_carrier_violations": len(baseline_violations),
        "baseline_carrier_max_ratio": baseline_violations[0]["ratio"],
        "top10_poisson_carrier_overlap": len(top10_nodes & violation_nodes),
        "top50_poisson_carrier_overlap": len(top50_nodes & violation_nodes),
        "poisson_carrier_log_correlation": pearson_log_correlation(
            baseline_residual_rows, baseline_ratios
        ),
        "minimum_candidate": min(candidate_rows[1:], key=lambda row: row["combined"]),
        "node_1257": next(
            row for row in top_poisson_with_steps if row["node_id"] == 1257
        ),
    }
    # Keep command stdout/stderr out of the compact result after retaining the
    # probe-specific numeric payload.
    for key in ("full_step", "block_steps"):
        summary[key].pop("stdout", None)
        summary[key].pop("stderr", None)
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
