#!/usr/bin/env python3
"""Recompute the 21-point TransportModels DD Id-Vg contact-basin baseline."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any

import transportmodels_fixed_contract as fixed


REPO = Path(__file__).resolve().parents[1]
UNCONFIGURED_ARTIFACT_ROOT = REPO / "__transportmodels_artifact_root_required__"
REF = UNCONFIGURED_ARTIFACT_ROOT
OLD_BASE = REF / "__unconfigured_old_base__.json"
OLD_SEED = REF / "__unconfigured_old_seed__.csv"
OLD_RUN = REF / "__unconfigured_old_run__"
OLD_ON_TAIL = REF / "__unconfigured_old_on_tail__"
SENTAURUS = REF / "__unconfigured_sentaurus__.csv"
OUTPUT = REF / "__unconfigured_output__"
REPORT_JSON = (
    REPO / "docs/validation/transportmodels_dd_contact_basin_v1_2026-08-24.json"
)
REPORT_MD = (
    REPO / "docs/validation/transportmodels_dd_contact_basin_v1_2026-08-24.md"
)
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"


def configure_artifact_paths(root: Path) -> None:
    """Bind all generated inputs and outputs to one explicit artifact bundle."""
    global REF, OLD_BASE, OLD_SEED, OLD_RUN, OLD_ON_TAIL, SENTAURUS, OUTPUT
    REF = root.resolve()
    OLD_BASE = (
        REF / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24/runs/dd"
        / "03_dd_idvg_curve.json"
    )
    OLD_SEED = OLD_BASE.with_name("dd_idvg_final_bias_relax_final_state.csv")
    OLD_RUN = OLD_BASE.parent
    OLD_ON_TAIL = OLD_BASE.parents[2] / "idvg_on_tail"
    SENTAURUS = REF / "run02/normalized/dd_idvg.csv"
    OUTPUT = REF / "vela_baseline/dd_contact_basin_fixed_contract_v1_2026-08-24"


def validate_artifact_bundle() -> None:
    required = [OLD_BASE, OLD_SEED, SENTAURUS]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "TransportModels artifact bundle is incomplete:\n  "
            + "\n  ".join(missing)
        )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def bias_tag(value: float) -> str:
    prefix = "m" if value < 0 else ""
    return prefix + f"{abs(value):.6f}".replace(".", "p")


def unique_bias_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_bias = {round(float(row["bias_V"]), 12): row for row in rows}
    return [by_bias[key] for key in sorted(by_bias)]


def environment() -> dict[str, str]:
    result = os.environ.copy()
    result["PATH"] = os.pathsep.join(
        [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin", result.get("PATH", "")]
    )
    return result


def exact_and_bridge_biases(exact: list[float]) -> list[float]:
    """Retain the 21-point report lattice and bridge weak inversion internally."""
    bridge_start = -0.52
    bridge_stop = -0.20
    count = int(round((bridge_stop - bridge_start) / 0.0025))
    bridge = [
        round(bridge_start + 0.0025 * index, 12)
        for index in range(1, count)
    ]
    return sorted(set(exact + bridge))


def execution_bounds(execution: list[float]) -> tuple[float, float] | None:
    """Return sweep bounds, or ``None`` when a resumed run is complete."""
    if not execution:
        return None
    return execution[0], execution[-1]


def configure() -> tuple[Path, list[float], list[float], float | None]:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    config = json.loads(OLD_BASE.read_text(encoding="utf-8"))
    config = fixed.apply_dd_contact_basin_contract(config)
    violations = fixed.validate_dd_contact_basin_config(config)
    if violations:
        raise RuntimeError(f"Contact-basin contract violations: {violations}")
    exact = fixed.load_contract()["bias_contract"]["idvg"]["gate_bias_V"]
    execution = exact_and_bridge_biases(exact)
    if any(right <= left for left, right in zip(execution, execution[1:])):
        raise RuntimeError("Contact-basin execution biases must be strictly increasing")
    raw_sources = [path for path in (
        OUTPUT / "dd_idvg_prefix.csv",
        OUTPUT / "dd_idvg_raw.csv",
        OUTPUT / "dd_idvg_segment.csv",
    ) if path.is_file()]
    raw_rows = unique_bias_rows([
        row for path in raw_sources for row in read_csv(path)
        if row.get("converged") == "1"
    ]) if raw_sources else []
    resume_bias = max((float(row["bias_V"]) for row in raw_rows), default=None)
    resume_state = None if resume_bias is None else OUTPUT / (
        f"state_bias_{bias_tag(resume_bias)}.csv"
    )
    if resume_state is not None and not resume_state.is_file():
        resume_bias = None
        resume_state = None
    resume = resume_bias is not None
    if resume:
        balance_sources = [path for path in (
            OUTPUT / "srh_balance_prefix.csv",
            OUTPUT / "srh_balance.csv",
            OUTPUT / "srh_balance_segment.csv",
        ) if path.is_file()]
        prefix_balance = unique_bias_rows([
            row for path in balance_sources for row in read_csv(path)
            if float(row["bias_V"]) <= resume_bias + 1e-12
        ])
        prefix_raw = [
            row for row in raw_rows if float(row["bias_V"]) <= resume_bias + 1e-12
        ]
        write_csv(OUTPUT / "dd_idvg_prefix.csv", prefix_raw)
        write_csv(OUTPUT / "srh_balance_prefix.csv", prefix_balance)
        execution = [bias for bias in execution if bias > resume_bias + 1e-12]
    path = OUTPUT / "config.json"
    bounds = execution_bounds(execution)
    if bounds is None:
        return path, exact, execution, resume_bias
    config["_comment"] = (
        "Frozen TransportModels DD contact-basin baseline; exact 21-point Id-Vg "
        "report lattice with an internal 2.5 mV weak-inversion bridge and "
        "block-filter globalization"
    )
    config["solver"].update({
        "line_search_mode": "block_filter",
        "residual_filter_gamma": 1.0e-4,
        "residual_filter_envelope_factor": 2.0,
        "continuity_row_scaling": {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": 1.0e-18,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12
        }
    })
    config["output_csv"] = str(
        (OUTPUT / ("dd_idvg_segment.csv" if resume else "dd_idvg_raw.csv")).resolve()
    )
    config["log_file"] = str((OUTPUT / "runner.log").resolve())
    sweep = config["sweep"]
    sweep.update({
        "start": bounds[0],
        "stop": bounds[1],
        "step": 0.0025,
        "bias_points": execution,
        "initial_state_file": str((resume_state if resume else OLD_SEED).resolve()),
        "write_state_file": str((OUTPUT / "final_state.csv").resolve()),
        "write_state_every_point_prefix": str((OUTPUT / "state").resolve()),
    })
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((OUTPUT / (
            "terminal_balance_segment.csv" if resume else "terminal_balance.csv"
        )).resolve()),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((OUTPUT / (
            "srh_balance_segment.csv" if resume else "srh_balance.csv"
        )).resolve()),
    }
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path, exact, execution, resume_bias


def merge_resumed_outputs() -> None:
    raw = read_csv(OUTPUT / "dd_idvg_prefix.csv")
    raw.extend(read_csv(OUTPUT / "dd_idvg_segment.csv"))
    write_csv(OUTPUT / "dd_idvg_raw.csv", raw)
    balance = read_csv(OUTPUT / "srh_balance_prefix.csv")
    balance.extend(read_csv(OUTPUT / "srh_balance_segment.csv"))
    write_csv(OUTPUT / "srh_balance.csv", balance)


def old_exact_state(bias: float) -> Path:
    current = OUTPUT / f"state_bias_{bias_tag(bias)}.csv"
    if current.is_file():
        return current
    if bias == -1.0:
        return OLD_SEED
    tag = bias_tag(bias)
    if bias <= 0.12:
        return OLD_RUN / f"dd_idvg_curve_state_bias_{tag}.csv"
    return OLD_ON_TAIL / f"state_bias_{tag}.csv"


def reclose_exact_points(
    runner: Path, exact: list[float]
) -> list[dict[str, Any]]:
    """Reclose every report point from the same-bias frozen physical state."""
    rows: list[dict[str, str]] = []
    balances: list[dict[str, str]] = []
    evidence: list[dict[str, Any]] = []
    sweep_curve_paths = [path for path in (
        OUTPUT / "dd_idvg_prefix.csv",
        OUTPUT / "dd_idvg_raw.csv",
        OUTPUT / "dd_idvg_segment.csv",
    ) if path.is_file()]
    sweep_balance_paths = [path for path in (
        OUTPUT / "srh_balance_prefix.csv",
        OUTPUT / "srh_balance.csv",
        OUTPUT / "srh_balance_segment.csv",
    ) if path.is_file()]
    sweep_rows = {
        round(float(row["bias_V"]), 12): row
        for path in sweep_curve_paths for row in read_csv(path)
        if row.get("converged") == "1"
    }
    sweep_balances = {
        round(float(row["bias_V"]), 12): row
        for path in sweep_balance_paths for row in read_csv(path)
    }
    for bias in exact:
        key = round(bias, 12)
        current_state = OUTPUT / f"state_bias_{bias_tag(bias)}.csv"
        if key in sweep_rows and key in sweep_balances and current_state.is_file():
            rows.append(sweep_rows[key])
            balances.append(sweep_balances[key])
            evidence.append({
                "bias_V": bias,
                "seed": "monotonic contact-basin sweep",
                "config": str((OUTPUT / "config.json").resolve()),
                "state": str(current_state.resolve()),
                "iterations": int(sweep_rows[key].get("iterations", "0")),
                "reused_completed_evidence": True,
            })
            continue
        seed = old_exact_state(bias)
        if not seed.is_file():
            raise FileNotFoundError(seed)
        run_dir = OUTPUT / "exact_reclosure" / bias_tag(bias)
        run_dir.mkdir(parents=True, exist_ok=True)
        existing_curve = run_dir / "curve.csv"
        existing_balance = run_dir / "srh_balance.csv"
        existing_state = run_dir / "final_state.csv"
        if existing_curve.is_file() and existing_balance.is_file() and existing_state.is_file():
            curve_rows = [
                row for row in read_csv(existing_curve) if row.get("converged") == "1"
            ]
            balance_rows = read_csv(existing_balance)
            if curve_rows and balance_rows:
                rows.append(curve_rows[-1])
                balances.append(balance_rows[-1])
                evidence.append({
                    "bias_V": bias,
                    "seed": "completed prior exact reclosure in this regression run",
                    "config": str((run_dir / "config.json").resolve()),
                    "state": str(existing_state.resolve()),
                    "iterations": int(curve_rows[-1].get("iterations", "0")),
                    "reused_completed_evidence": True,
                })
                continue
        config = json.loads(OLD_BASE.read_text(encoding="utf-8"))
        config = fixed.apply_dd_contact_basin_contract(config)
        for contact in config["contacts"]:
            if contact["name"] == "gate":
                contact["bias"] = bias
        config["solver"].update({
            "line_search_mode": "block_filter",
            "residual_filter_gamma": 1.0e-4,
            "residual_filter_envelope_factor": 2.0,
            "continuity_row_scaling": {
                "enabled": True,
                "flux_fraction": 1.0e-3,
                "scale_floor": 1.0e-30,
                "min_source_scale": 1.0e-18,
                "min_weight": 1.0e-12,
                "max_weight": 1.0e12
            }
        })
        config["output_csv"] = str((run_dir / "curve.csv").resolve())
        config["log_file"] = str((run_dir / "runner.log").resolve())
        sweep = config["sweep"]
        sweep.pop("initialization", None)
        sweep.update({
            "start": bias,
            "stop": bias,
            "step": 0.01,
            "bias_points": [bias],
            "initial_state_file": str(seed.resolve()),
            "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        })
        diagnostics = sweep.setdefault("diagnostics", {})
        diagnostics["terminal_balance"] = {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str((run_dir / "terminal_balance.csv").resolve()),
        }
        diagnostics["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str((run_dir / "srh_balance.csv").resolve()),
        }
        config_path = run_dir / "config.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)], cwd=REPO,
            env=environment(), text=True, capture_output=True, check=False,
        )
        (run_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (run_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
        curve_rows = [
            row for row in read_csv(run_dir / "curve.csv")
            if row.get("converged") == "1"
        ]
        if completed.returncode or not curve_rows:
            raise RuntimeError(f"Contact-basin exact reclosure failed at Vg={bias}")
        rows.append(curve_rows[-1])
        balance_rows = read_csv(run_dir / "srh_balance.csv")
        if not balance_rows:
            raise RuntimeError(f"Missing SRH/KCL evidence at Vg={bias}")
        balances.append(balance_rows[-1])
        evidence.append({
            "bias_V": bias,
            "seed": str(seed.resolve()),
            "config": str(config_path.resolve()),
            "state": str((run_dir / "final_state.csv").resolve()),
            "iterations": int(curve_rows[-1].get("iterations", "0")),
        })
    write_csv(OUTPUT / "dd_idvg_raw.csv", rows)
    write_csv(OUTPUT / "srh_balance.csv", balances)
    return evidence


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def align(exact: list[float]) -> list[dict[str, Any]]:
    raw = {
        round(float(row["bias_V"]), 12): row
        for row in read_csv(OUTPUT / "dd_idvg_raw.csv")
        if row.get("converged") == "1"
    }
    balance = {
        round(float(row["bias_V"]), 12): row
        for row in read_csv(OUTPUT / "srh_balance.csv")
    }
    reference = {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in read_csv(SENTAURUS)
    }
    missing = [bias for bias in exact if round(bias, 12) not in raw]
    if missing:
        raise RuntimeError(f"Missing exact DD bias points: {missing}")
    candidate_rows: list[dict[str, Any]] = []
    aligned: list[dict[str, Any]] = []
    for bias in exact:
        key = round(bias, 12)
        row = raw[key]
        candidate_rows.append(row)
        vela = abs(float(row["current_total_A_per_um"]))
        sentaurus = reference[key]
        balance_row = balance.get(key, {})
        aligned.append({
            "bias_V": bias,
            "vela_A_per_um": vela,
            "sentaurus_A_per_um": sentaurus,
            "absolute_relative_error": abs(vela - sentaurus) / max(sentaurus, 1e-300),
            "absolute_log_error_dex": abs(
                math.log10(max(vela, 1e-300))
                - math.log10(max(sentaurus, 1e-300))
            ),
            "four_terminal_kcl_residual_A_per_um": abs(float(
                balance_row.get("four_terminal_kcl_residual_A_per_um", "nan")
            )),
            "id_to_kcl_residual_ratio": float(
                balance_row.get("id_to_kcl_residual_ratio", "nan")
            ),
            "numerical_status": balance_row.get("numerical_status", "missing"),
        })
    write_csv(OUTPUT / "dd_idvg_21_point.csv", candidate_rows)
    write_csv(OUTPUT / "dd_idvg_21_point_aligned.csv", aligned)
    return aligned


def summarize(aligned: list[dict[str, Any]]) -> dict[str, Any]:
    limits = fixed.load_dd_contact_basin_contract()["acceptance"]
    regions = {
        "deep_off": aligned[:3],
        "transition": aligned[3:8],
        "on": aligned[8:],
    }
    metrics = {
        name: {
            "max_absolute_relative_error": max(
                row["absolute_relative_error"] for row in rows
            ),
            "max_absolute_log_error_dex": max(
                row["absolute_log_error_dex"] for row in rows
            ),
            "median_absolute_log_error_dex": percentile(
                [row["absolute_log_error_dex"] for row in rows], 0.5
            ),
        }
        for name, rows in regions.items()
    }
    deep_points = []
    for row in regions["deep_off"]:
        log_pass = (
            row["absolute_log_error_dex"]
            <= limits["deep_off_max_absolute_log_error_dex"]
        )
        resolved = (
            row["id_to_kcl_residual_ratio"]
            >= limits["deep_off_min_id_to_kcl_ratio"]
            and row["numerical_status"] != "numerically_unresolved"
        )
        deep_points.append({
            **row,
            "status": "pass" if log_pass and resolved else
                      "numerically_unresolved" if not resolved else "fail",
        })
    main_pass = (
        metrics["transition"]["max_absolute_log_error_dex"]
        <= limits["idvg_transition_max_absolute_log_error_dex"]
        and metrics["on"]["max_absolute_relative_error"]
        <= limits["idvg_on_max_absolute_relative_error"]
    )
    deep_pass = all(row["status"] == "pass" for row in deep_points)
    return {
        "metrics": metrics,
        "acceptance": {
            "main_curve_pass": main_pass,
            "deep_off_pass": deep_pass,
            "overall_pass": main_pass and deep_pass,
            "deep_off_points": deep_points,
        },
    }


def make_plot(aligned: list[dict[str, Any]]) -> str | None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    bias = [row["bias_V"] for row in aligned]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    axes[0].semilogy(bias, [row["sentaurus_A_per_um"] for row in aligned],
                    "-", label="Sentaurus 2022")
    axes[0].semilogy(bias, [row["vela_A_per_um"] for row in aligned],
                    "o", ms=3.5, label="Vela DD contact_basin")
    axes[0].set(xlabel="Gate voltage Vg (V)", ylabel="Drain current Id (A/um)",
                title="21-point DD Id-Vg")
    axes[1].plot(bias, [100 * row["absolute_relative_error"] for row in aligned],
                 "o-", ms=3.5)
    axes[1].set(xlabel="Gate voltage Vg (V)", ylabel="Absolute relative error (%)",
                title="Pointwise relative error")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
    axes[0].legend()
    fig.tight_layout()
    path = OUTPUT / "dd_idvg_21_point_comparison.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path.resolve())


def write_report(report: dict[str, Any]) -> None:
    payload, nonfinite_paths = fixed.strict_json_payload(report)
    if nonfinite_paths:
        payload["serialization"] = {
            "nonfinite_values_replaced_with_null": nonfinite_paths,
            "reason": "source evidence did not provide a finite value",
        }
    REPORT_JSON.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    deep_rows = []
    for row in report["acceptance"]["deep_off_points"]:
        deep_rows.append(
            f"| {row['bias_V']:.2f} | {row['sentaurus_A_per_um']:.6e} | "
            f"{row['vela_A_per_um']:.6e} | {100*row['absolute_relative_error']:.3f}% | "
            f"{row['absolute_log_error_dex']:.6f} | "
            f"{row['id_to_kcl_residual_ratio']:.3f} | {row['status']} |"
        )
    metrics = report["metrics"]
    markdown = f"""# TransportModels DD contact_basin 21点 Id-Vg 回归

## 结论

- 计算状态：完成，共21个固定比较点，另使用 {report['execution']['internal_bridge_points']} 个内部延续点。
- 主曲线验收：{'通过' if report['acceptance']['main_curve_pass'] else '未通过'}。
- 深关断数值验收：{'通过' if report['acceptance']['deep_off_pass'] else '未通过'}。
- 总体验收：{'通过' if report['acceptance']['overall_pass'] else '未通过'}。

## 分区误差

| 区域 | 最大相对误差 | 最大对数误差 |
|---|---:|---:|
| 深关断 | {100*metrics['deep_off']['max_absolute_relative_error']:.3f}% | {metrics['deep_off']['max_absolute_log_error_dex']:.6f} dex |
| 过渡区 | {100*metrics['transition']['max_absolute_relative_error']:.3f}% | {metrics['transition']['max_absolute_log_error_dex']:.6f} dex |
| 导通区 | {100*metrics['on']['max_absolute_relative_error']:.3f}% | {metrics['on']['max_absolute_log_error_dex']:.6f} dex |

## 深关断前三点

| Vg (V) | Sentaurus Id (A/um) | Vela Id (A/um) | 相对误差 | 对数误差 (dex) | Id/abs(KCL) | 状态 |
|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(deep_rows)}

## 固定证据

- 物理契约：`{report['contracts']['base']['path']}`
- DD 数值契约：`{report['contracts']['contact_basin']['path']}`
- 运行配置：`{report['execution']['config']}`
- 21点曲线：`{report['artifacts']['curve_csv']}`
- 对齐结果：`{report['artifacts']['aligned_csv']}`
"""
    REPORT_MD.write_text(markdown, encoding="utf-8")


def acceptance_exit_code(acceptance: dict[str, Any]) -> int:
    """Make regression failure observable to CTest and other automation."""
    return 0 if acceptance.get("overall_pass") is True else 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=None,
        help=(
            "Mounted TransportModels artifact bundle. Alternatively set "
            f"{fixed.ARTIFACT_ROOT_ENV}."
        ),
    )
    args = parser.parse_args()
    try:
        configure_artifact_paths(fixed.resolve_artifact_root(args.artifact_root))
        validate_artifact_bundle()
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(runner)
    if not OLD_BASE.is_file() or not OLD_SEED.is_file():
        raise FileNotFoundError("The frozen DD base config or -1 V seed is missing")
    attempts: list[dict[str, Any]] = []
    exact: list[float] = []
    execution: list[float] = []
    resume_bias: float | None = None
    resumed = False
    config_path = OUTPUT / "config.json"
    sweep_completed = False
    already_complete = False
    stalled = False
    for attempt in range(1, 21):
        config_path, exact, execution, resume_bias = configure()
        resumed = resumed or resume_bias is not None
        previous_bias = -math.inf if resume_bias is None else resume_bias
        if not execution:
            sweep_completed = True
            already_complete = True
            break
        (OUTPUT / f"config_attempt_{attempt:02d}.json").write_text(
            config_path.read_text(encoding="utf-8"), encoding="utf-8"
        )
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)],
            cwd=REPO, env=environment(), text=True, capture_output=True, check=False,
        )
        (OUTPUT / f"stdout_attempt_{attempt:02d}.txt").write_text(
            completed.stdout, encoding="utf-8"
        )
        (OUTPUT / f"stderr_attempt_{attempt:02d}.txt").write_text(
            completed.stderr, encoding="utf-8"
        )
        segment = [
            row for row in read_csv(OUTPUT / (
                "dd_idvg_segment.csv" if resume_bias is not None else "dd_idvg_raw.csv"
            )) if row.get("converged") == "1"
        ]
        latest_bias = max(
            (float(row["bias_V"]) for row in segment), default=previous_bias
        )
        attempts.append({
            "attempt": attempt,
            "resume_bias_V": None if not math.isfinite(previous_bias) else previous_bias,
            "last_converged_bias_V": latest_bias,
            "returncode": completed.returncode,
            "made_progress": latest_bias > previous_bias + 1e-12,
        })
        if completed.returncode == 0:
            sweep_completed = True
            break
        if latest_bias <= previous_bias + 1e-12:
            stalled = True
            break
    else:
        stalled = True
    reclosure_evidence: list[dict[str, Any]] = []
    execution_mode = "monotonic_sweep"
    if sweep_completed and resumed and not already_complete:
        merge_resumed_outputs()
    elif not sweep_completed:
        execution_mode = "same_bias_exact_reclosure_after_sweep_stall"
        reclosure_evidence = reclose_exact_points(runner, exact)
    aligned = align(exact)
    summary = summarize(aligned)
    plot = make_plot(aligned)
    numerical_contract = fixed.DEFAULT_DD_CONTACT_BASIN_CONTRACT
    base_contract = fixed.DEFAULT_CONTRACT
    report = {
        "schema": "vela.transportmodels.dd.contact_basin.acceptance.v1",
        "as_of": "2026-08-24",
        "execution": {
            "status": "complete",
            "mode": execution_mode,
            "reported_points": len(exact),
            "executed_points": len(exact_and_bridge_biases(exact)),
            "current_process_points": len(execution),
            "resumed": resumed,
            "resume_bias_V": resume_bias,
            "progress_restart_attempts": attempts,
            "sweep_stalled": stalled,
            "exact_reclosures": reclosure_evidence,
            "internal_bridge_points": len(exact_and_bridge_biases(exact)) - len(exact),
            "config": str(config_path.resolve()),
            "runner": str(runner),
            "runner_sha256": fixed.sha256(runner),
        },
        "contracts": {
            "base": {"path": str(base_contract.resolve()),
                     "sha256": fixed.sha256(base_contract)},
            "contact_basin": {"path": str(numerical_contract.resolve()),
                              "sha256": fixed.sha256(numerical_contract)},
        },
        "aligned": aligned,
        **summary,
        "artifacts": {
            "curve_csv": str((OUTPUT / "dd_idvg_21_point.csv").resolve()),
            "aligned_csv": str((OUTPUT / "dd_idvg_21_point_aligned.csv").resolve()),
            "plot_png": plot,
            "raw_curve_csv": str((OUTPUT / "dd_idvg_raw.csv").resolve()),
            "srh_balance_csv": str((OUTPUT / "srh_balance.csv").resolve()),
        },
    }
    write_report(report)
    terminal_payload = {
        "status": "complete",
        "reported_points": len(exact),
        "metrics": summary["metrics"],
        "acceptance": summary["acceptance"],
        "report": str(REPORT_JSON.resolve()),
    }
    print(fixed.strict_json_text(terminal_payload))
    return acceptance_exit_code(summary["acceptance"])


if __name__ == "__main__":
    raise SystemExit(main())
