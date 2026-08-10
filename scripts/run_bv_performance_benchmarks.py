#!/usr/bin/env python3
"""Prepare and run the three controlled BV performance benchmarks.

The benchmarks deliberately reuse only checkpoints that precede the measured
work.  Each run gets an isolated directory so boundary-control persistence and
gmon.out cannot leak between repetitions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import time
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
VALIDATION = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation"
)
EXTERNAL_SOURCE = VALIDATION / "boundary_external_resistor_20260806"
CURRENT_SOURCE = VALIDATION / "boundary_voltage_to_current_20260806"
HIGH_FIELD_SOURCE = (
    EXTERNAL_SOURCE / "repro_6p087_to_6p099_diagfix_max220_20260807"
)

SCENARIOS = (
    "high_field_transition",
    "voltage_to_current_final",
    "external_resistor_1206",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def command_output(command: list[str]) -> str:
    return subprocess.check_output(command, text=True, encoding="utf-8").strip()


def physical_memory_bytes() -> int | None:
    if os.name != "nt":
        try:
            return int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
        except (AttributeError, OSError, ValueError):
            return None
    try:
        import ctypes

        class MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("length", ctypes.c_ulong),
                ("memory_load", ctypes.c_ulong),
                ("total_physical", ctypes.c_ulonglong),
                ("available_physical", ctypes.c_ulonglong),
                ("total_page_file", ctypes.c_ulonglong),
                ("available_page_file", ctypes.c_ulonglong),
                ("total_virtual", ctypes.c_ulonglong),
                ("available_virtual", ctypes.c_ulonglong),
                ("available_extended_virtual", ctypes.c_ulonglong),
            ]

        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return int(status.total_physical)
    except (AttributeError, OSError):
        pass
    return None


def cpu_name() -> str:
    if os.name == "nt":
        try:
            value = command_output(
                [
                    "powershell.exe",
                    "-NoProfile",
                    "-Command",
                    "(Get-ItemProperty "
                    "'HKLM:\\HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0'"
                    ").ProcessorNameString",
                ]
            )
            if value:
                return value
        except (OSError, subprocess.CalledProcessError):
            pass
    return platform.processor() or "unknown"


def configure_output_paths(
    config: dict[str, Any],
    output: Path,
    enable_performance_profiling: bool = False,
) -> None:
    config["output_csv"] = str((output / "sweep.csv").resolve())
    solver = config["solver"]
    if enable_performance_profiling:
        solver["performance_profiling"] = {
            "enabled": True,
            "json_file": str((output / "performance_profile.json").resolve()),
        }
    carrier = solver.get("carrier_row_convergence", {})
    if "diagnostic_csv" in carrier:
        carrier["diagnostic_csv"] = str(
            (output / "carrier_row_convergence.csv").resolve()
        )
    if "trace_csv" in carrier:
        carrier["trace_csv"] = str((output / "carrier_row_trace.csv").resolve())

    sweep = config["sweep"]
    sweep["write_state_file"] = str((output / "last_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (output / "states" / "accepted_state").resolve()
    )
    diagnostics = sweep.get("diagnostics", {})
    if "qf_bounds" in diagnostics:
        diagnostics["qf_bounds"]["csv_file"] = str(
            (output / "qf_bounds.csv").resolve()
        )
    history = diagnostics.get("newton_history")
    if history:
        history["csv_file"] = str((output / "newton_history.csv").resolve())
        history["attempts_csv_file"] = str(
            (output / "newton_attempts.csv").resolve()
        )
        history["iterations_csv_file"] = str(
            (output / "newton_iterations.csv").resolve()
        )
    control = sweep.get("boundary_control")
    if control:
        control["evaluation_csv"] = str(
            (output / "boundary_control_evaluations.csv").resolve()
        )
        control["checkpoint_directory"] = str(
            (output / "boundary_control_checkpoints").resolve()
        )


def copy_initial_state(config: dict[str, Any], output: Path) -> None:
    source = Path(config["sweep"]["initial_state_file"])
    destination = output / "inputs" / "initial_state.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    config["sweep"]["initial_state_file"] = str(destination.resolve())


def selected_seed_rows(
    scenario: str, source_csv: Path
) -> tuple[list[str], list[dict[str, str]]]:
    with source_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = []
        for row in reader:
            target = float(row["target_value"])
            evaluation = int(row["evaluation_index"])
            if keep_seed_evaluation(scenario, target, evaluation):
                rows.append(row)
    return fieldnames, rows


def keep_seed_evaluation(scenario: str, target: float, evaluation: int) -> bool:
    if scenario == "voltage_to_current_final":
        return target < 1.0e-4 or (
            abs(target - 1.0e-4) <= 1.0e-15 and evaluation <= 6
        )
    if scenario == "external_resistor_1206":
        return target < 1206.0
    return False


def copy_boundary_seed(
    scenario: str, source: Path, output: Path
) -> int:
    fieldnames, rows = selected_seed_rows(
        scenario, source / "boundary_control_evaluations.csv"
    )
    checkpoints = output / "boundary_control_checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    copied: dict[Path, Path] = {}
    for row in rows:
        state_text = row.get("state_file", "")
        if not state_text:
            continue
        state_source = Path(state_text)
        if not state_source.exists():
            raise FileNotFoundError(f"missing seed checkpoint: {state_source}")
        state_destination = checkpoints / state_source.name
        if state_source not in copied:
            shutil.copy2(state_source, state_destination)
            copied[state_source] = state_destination
        row["state_file"] = str(state_destination.resolve())

    evaluation_csv = output / "boundary_control_evaluations.csv"
    with evaluation_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def scenario_source(scenario: str) -> Path:
    if scenario == "high_field_transition":
        return HIGH_FIELD_SOURCE
    if scenario == "voltage_to_current_final":
        return CURRENT_SOURCE
    if scenario == "external_resistor_1206":
        return EXTERNAL_SOURCE
    raise ValueError(f"unknown scenario: {scenario}")


def configure_scenario_optimizations(
    scenario: str, config: dict[str, Any]
) -> None:
    if scenario != "external_resistor_1206":
        return
    config["sweep"]["boundary_control"]["predictor_max_step_factor"] = 3.0
    config["sweep"]["continuation"] = {
        "predictor": {
            "mode": "secant",
            "fields": ["psi", "phin", "phip"],
            "max_extrapolation_ratio": 4.0,
        }
    }


def prepare_scenario(
    scenario: str,
    output: Path,
    enable_performance_profiling: bool = False,
) -> dict[str, Any]:
    source = scenario_source(scenario)
    config_source = source / "simulation.json"
    config = load_json(config_source)
    output.mkdir(parents=True, exist_ok=True)
    configure_output_paths(config, output, enable_performance_profiling)
    copy_initial_state(config, output)

    seed_rows = 0
    if scenario == "high_field_transition":
        config["sweep"]["boundary_control"]["resume"] = False
    else:
        seed_rows = copy_boundary_seed(scenario, source, output)
    configure_scenario_optimizations(scenario, config)

    config_path = output / "simulation.json"
    write_json(config_path, config)
    physics = {
        "mobility": config["solver"].get("mobility"),
        "impact_ionization": config["solver"].get("impact_ionization"),
        "band_to_band": config["solver"].get("band_to_band"),
        "carrier_statistics": config["solver"].get("carrier_statistics"),
    }
    manifest = {
        "schema_version": 1,
        "scenario": scenario,
        "source_config": str(config_source.resolve()),
        "source_config_sha256": sha256_file(config_source),
        "prepared_config_sha256": sha256_file(config_path),
        "physics_config_sha256": sha256_json(physics),
        "initial_state_sha256": sha256_file(output / "inputs/initial_state.csv"),
        "seed_evaluation_rows": seed_rows,
    }
    write_json(output / "benchmark_manifest.json", manifest)
    return manifest


def read_appended_boundary_rows(output: Path, seed_rows: int) -> list[dict[str, str]]:
    path = output / "boundary_control_evaluations.csv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return rows[seed_rows:]


def last_csv_row(path: Path) -> dict[str, str] | None:
    if not path.exists():
        return None
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return rows[-1] if rows else None


def extract_bv(
    sweep_csv: Path,
    threshold: float = 1.0e-4,
    current_tolerance: float = 0.0,
) -> float | None:
    if not sweep_csv.exists():
        return None
    with sweep_csv.open(newline="", encoding="utf-8-sig") as handle:
        rows = [row for row in csv.DictReader(handle) if row["converged"] == "1"]
    pairs = [
        (float(row["inner_voltage_V"]), abs(float(row["current_total_A_per_um"])))
        for row in rows
        if row.get("inner_voltage_V") and row.get("current_total_A_per_um")
    ]
    for (voltage0, current0), (voltage1, current1) in zip(pairs, pairs[1:]):
        if current0 <= threshold <= current1:
            if current1 == current0:
                return voltage1
            fraction = (threshold - current0) / (current1 - current0)
            return voltage0 + fraction * (voltage1 - voltage0)
    if current_tolerance > 0.0 and pairs:
        voltage, current = min(
            pairs, key=lambda pair: abs(pair[1] - threshold)
        )
        if abs(current - threshold) <= current_tolerance:
            return voltage
    return None


def parse_result(output: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    scenario = manifest["scenario"]
    if scenario == "high_field_transition":
        row = last_csv_row(output / "sweep.csv")
        return {
            "full_dd_evaluations": 1 if row else 0,
            "newton_iterations": int(row["newton_iterations"]) if row else 0,
            "converged": row["converged"] == "1" if row else False,
            "inner_voltage_V": float(row["bias_V"]) if row else None,
            "global_electron_continuity_ratio": (
                float(row["global_electron_continuity_closure_ratio"])
                if row and row["global_electron_continuity_closure_ratio"]
                else None
            ),
            "global_hole_continuity_ratio": (
                float(row["global_hole_continuity_closure_ratio"])
                if row and row["global_hole_continuity_closure_ratio"]
                else None
            ),
        }

    appended = read_appended_boundary_rows(
        output, int(manifest["seed_evaluation_rows"])
    )
    fresh = [row for row in appended if row.get("resumed") == "0"]
    sweep_row = last_csv_row(output / "sweep.csv")
    config = json.loads((output / "simulation.json").read_text(encoding="utf-8"))
    validation = config.get("_validation_case", {})
    voltage_to_current = config.get("sweep", {}).get("voltage_to_current", {})
    bv_threshold = float(validation.get("current_threshold_A_per_um", 1.0e-4))
    current_tolerance = float(
        voltage_to_current.get("current_tolerance_A_per_um", 0.0)
    )
    result: dict[str, Any] = {
        "full_dd_evaluations": len(fresh),
        "newton_iterations": sum(int(row["newton_iterations"]) for row in fresh),
        "converged": bool(sweep_row and sweep_row["converged"] == "1"),
        "inner_voltage_V": (
            float(sweep_row["inner_voltage_V"])
            if sweep_row and sweep_row["inner_voltage_V"]
            else None
        ),
        "global_electron_continuity_ratio": (
            float(sweep_row["global_electron_continuity_closure_ratio"])
            if sweep_row and sweep_row["global_electron_continuity_closure_ratio"]
            else None
        ),
        "global_hole_continuity_ratio": (
            float(sweep_row["global_hole_continuity_closure_ratio"])
            if sweep_row and sweep_row["global_hole_continuity_closure_ratio"]
            else None
        ),
        "vela_bv_V": extract_bv(
            output / "sweep.csv", bv_threshold, current_tolerance
        ),
    }
    if sweep_row:
        if scenario == "voltage_to_current_final":
            result["boundary_residual"] = float(
                sweep_row["current_boundary_residual_A_per_um"]
            )
            result["boundary_residual_unit"] = "A_per_um"
        else:
            result["boundary_residual"] = float(sweep_row["load_line_residual_V"])
            result["boundary_residual_unit"] = "V"
    return result


def build_metadata(executable: Path, command: list[str]) -> dict[str, Any]:
    git = command_output(["git", "-C", str(REPO), "rev-parse", "HEAD"])
    return {
        "git_commit": git,
        "executable": str(executable.resolve()),
        "executable_sha256": sha256_file(executable),
        "command": command,
        "environment": {
            "cpu": cpu_name(),
            "logical_processors": os.cpu_count(),
            "physical_memory_bytes": physical_memory_bytes(),
            "os": platform.platform(),
            "python": platform.python_version(),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
    }


def parse_gprof_flat(flat_path: Path, csv_path: Path) -> int:
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    with_calls = re.compile(
        rf"^\s*({number})\s+({number})\s+({number})\s+(\d+)\s+"
        rf"({number})\s+({number})\s+(.+?)\s*$"
    )
    without_calls = re.compile(
        rf"^\s*({number})\s+({number})\s+({number})\s+(.+?)\s*$"
    )
    rows: list[dict[str, str]] = []
    for line in flat_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = with_calls.match(line)
        if match:
            percent, cumulative, self_seconds, calls, self_call, total_call, name = (
                match.groups()
            )
        else:
            match = without_calls.match(line)
            if not match:
                continue
            percent, cumulative, self_seconds, name = match.groups()
            calls = self_call = total_call = ""
        if name.startswith("time ") or name.startswith("name "):
            continue
        rows.append(
            {
                "percent_time": percent,
                "cumulative_seconds": cumulative,
                "self_seconds": self_seconds,
                "calls": calls,
                "self_seconds_per_call": self_call,
                "total_seconds_per_call": total_call,
                "function": name,
            }
        )
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "percent_time",
                "cumulative_seconds",
                "self_seconds",
                "calls",
                "self_seconds_per_call",
                "total_seconds_per_call",
                "function",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def parse_gprof_callgraph(callgraph_path: Path, csv_path: Path) -> int:
    """Extract the primary, cumulative-time-ranked gprof call-graph rows."""
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    primary = re.compile(
        rf"^\[(\d+)\]\s+({number})\s+({number})\s+({number})\s+"
        rf"(?:(\d+)\s+)?(.+?)\s+\[\d+\]\s*$"
    )
    rows: list[dict[str, str]] = []
    for line in callgraph_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        match = primary.match(line)
        if not match:
            continue
        index, percent, self_seconds, children_seconds, calls, name = match.groups()
        rows.append(
            {
                "index": index,
                "percent_time": percent,
                "self_seconds": self_seconds,
                "children_seconds": children_seconds,
                "cumulative_seconds": str(
                    float(self_seconds) + float(children_seconds)
                ),
                "calls": calls or "",
                "function": name,
            }
        )
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "index",
                "percent_time",
                "self_seconds",
                "children_seconds",
                "cumulative_seconds",
                "calls",
                "function",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def summarize_gprof(output: Path) -> dict[str, Any] | None:
    flat_path = output / "gprof_hotspots.csv"
    callgraph_path = output / "gprof_callgraph_hotspots.csv"
    if not flat_path.exists() or not callgraph_path.exists():
        return None

    def read_rows(path: Path) -> list[dict[str, str]]:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))

    flat = read_rows(flat_path)
    callgraph = read_rows(callgraph_path)
    profiler_runtime = {"_mcount_private", "__fentry__"}
    production_flat = [row for row in flat if row["function"] not in profiler_runtime]
    production_callgraph = [
        row for row in callgraph if row["function"] not in profiler_runtime
    ]
    calls_ranked = sorted(
        (row for row in production_flat if row.get("calls")),
        key=lambda row: int(row["calls"]),
        reverse=True,
    )
    artifacts: dict[str, Any] = {}
    for name in (
        "gmon.out",
        "gprof_flat.txt",
        "gprof_callgraph.txt",
        "gprof_hotspots.csv",
        "gprof_callgraph_hotspots.csv",
    ):
        path = output / name
        if path.exists():
            artifacts[name] = {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    return {
        "artifacts": artifacts,
        "profiler_runtime_percent": sum(
            float(row["percent_time"])
            for row in flat
            if row["function"] in profiler_runtime
        ),
        "top_self_time": production_flat[:15],
        "top_cumulative_time": production_callgraph[:15],
        "top_call_counts": calls_ranked[:15],
        "self_time_candidates_over_5_percent": [
            row for row in production_flat if float(row["percent_time"]) >= 5.0
        ],
        "cumulative_time_candidates_over_10_percent": [
            row
            for row in production_callgraph
            if float(row["percent_time"]) >= 10.0
        ],
    }


def generate_gprof(output: Path, executable: Path, gprof: Path) -> None:
    gmon = output / "gmon.out"
    if not gmon.exists():
        raise FileNotFoundError(f"profile run did not produce {gmon}")
    flat = output / "gprof_flat.txt"
    callgraph = output / "gprof_callgraph.txt"
    with flat.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [str(gprof), "-b", "-p", str(executable), str(gmon)],
            check=True,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    with callgraph.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [str(gprof), "-b", "-q", str(executable), str(gmon)],
            check=True,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    hotspot_rows = parse_gprof_flat(flat, output / "gprof_hotspots.csv")
    if hotspot_rows == 0:
        raise RuntimeError(
            "gprof produced no mapped hotspot rows; check profiling linkage and ASLR"
        )
    callgraph_rows = parse_gprof_callgraph(
        callgraph, output / "gprof_callgraph_hotspots.csv"
    )
    if callgraph_rows == 0:
        raise RuntimeError("gprof produced no mapped call-graph rows")


def run_scenario(
    scenario: str,
    output: Path,
    executable: Path,
    gprof: Path | None,
    runtime_bin: Path | None,
    enable_performance_profiling: bool = False,
) -> dict[str, Any]:
    manifest = prepare_scenario(
        scenario, output, enable_performance_profiling
    )
    command = [str(executable.resolve()), "--config", str((output / "simulation.json").resolve())]
    metadata = build_metadata(executable, command)
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    if runtime_bin is not None:
        environment["PATH"] = str(runtime_bin.resolve()) + os.pathsep + environment.get(
            "PATH", ""
        )
        performance_runtime_bin = str(runtime_bin.resolve())
    else:
        performance_runtime_bin = None
    start_wall = time.time()
    start = time.perf_counter()
    with (output / "run.log").open("w", encoding="utf-8") as log:
        process = subprocess.run(
            command,
            cwd=output,
            env=environment,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    elapsed = time.perf_counter() - start
    performance = {
        "schema_version": 1,
        "scenario": scenario,
        "started_unix_s": start_wall,
        "wall_seconds": elapsed,
        "exit_code": process.returncode,
        **metadata,
        "runtime_bin": performance_runtime_bin,
        "benchmark": manifest,
        "result": parse_result(output, manifest),
    }
    write_json(output / "performance_run.json", performance)
    if gprof is not None and process.returncode == 0:
        generate_gprof(output, executable, gprof)
    return performance


def aggregate_runs(root: Path, scenarios: Iterable[str]) -> dict[str, Any]:
    aggregate: dict[str, Any] = {"schema_version": 1, "scenarios": {}}
    phase_rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        runs = []
        for path in sorted((root / scenario).glob("rep_*/performance_run.json")):
            run = load_json(path)
            run["result"] = parse_result(path.parent, run["benchmark"])
            profile_path = path.parent / "performance_profile.json"
            if profile_path.exists():
                run["internal_profile"] = load_json(profile_path)
            write_json(path, run)
            runs.append(run)
        times = [float(run["wall_seconds"]) for run in runs if run["exit_code"] == 0]
        median = statistics.median(times) if times else None
        max_deviation = (
            max(abs(value - median) for value in times) if median is not None else None
        )
        aggregate["scenarios"][scenario] = {
            "successful_runs": len(times),
            "wall_seconds": times,
            "median_wall_seconds": median,
            "max_deviation_seconds": max_deviation,
            "max_deviation_fraction": (
                max_deviation / median if median and max_deviation is not None else None
            ),
            "results": [run["result"] for run in runs],
            "internal_profiles": [
                run["internal_profile"]
                for run in runs
                if "internal_profile" in run
            ],
        }
        if runs:
            gprof = summarize_gprof(
                root / scenario / "rep_01"
            )
            if gprof is not None:
                aggregate["scenarios"][scenario]["gprof"] = gprof
        for index, run in enumerate(runs, start=1):
            phase_rows.append(
                {
                    "scenario": scenario,
                    "repetition": index,
                    "wall_seconds": run["wall_seconds"],
                    "exit_code": run["exit_code"],
                    "full_dd_evaluations": run["result"]["full_dd_evaluations"],
                    "newton_iterations": run["result"]["newton_iterations"],
                }
            )
    write_json(root / "performance_summary.json", aggregate)
    with (root / "performance_phases.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scenario",
                "repetition",
                "wall_seconds",
                "exit_code",
                "full_dd_evaluations",
                "newton_iterations",
            ],
        )
        writer.writeheader()
        writer.writerows(phase_rows)
    return aggregate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", choices=(*SCENARIOS, "all"), default="all")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--gprof", type=Path)
    parser.add_argument(
        "--enable-performance-profiling",
        action="store_true",
        help="enable low-overhead internal stage timers and counters",
    )
    parser.add_argument(
        "--runtime-bin",
        type=Path,
        default=Path(r"D:\msys64\ucrt64\bin") if os.name == "nt" else None,
        help="directory containing compiler runtime DLLs",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="write isolated inputs and manifests without running the solver",
    )
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        help="refresh result parsing and aggregate existing run directories",
    )
    args = parser.parse_args()

    scenarios = SCENARIOS if args.scenario == "all" else (args.scenario,)
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.summarize_only:
        print(json.dumps(aggregate_runs(args.output_root, SCENARIOS), indent=2))
        return 0
    if args.runner is None:
        parser.error("--runner is required unless --summarize-only is used")
    failures = 0
    for scenario in scenarios:
        for repetition in range(1, args.repetitions + 1):
            output = args.output_root / scenario / f"rep_{repetition:02d}"
            if output.exists() and any(output.iterdir()):
                raise FileExistsError(
                    f"refusing to reuse non-empty benchmark directory: {output}"
                )
            print(f"running {scenario} repetition {repetition}: {output}", flush=True)
            if args.prepare_only:
                manifest = prepare_scenario(
                    scenario, output, args.enable_performance_profiling
                )
                print(json.dumps(manifest, indent=2), flush=True)
                continue
            performance = run_scenario(
                scenario,
                output,
                args.runner.resolve(),
                args.gprof,
                args.runtime_bin,
                args.enable_performance_profiling,
            )
            print(
                json.dumps(
                    {
                        "scenario": scenario,
                        "repetition": repetition,
                        "wall_seconds": performance["wall_seconds"],
                        "exit_code": performance["exit_code"],
                        "result": performance["result"],
                    },
                    indent=2,
                ),
                flush=True,
            )
            failures += performance["exit_code"] != 0
    if not args.prepare_only:
        aggregate_runs(args.output_root, SCENARIOS)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
