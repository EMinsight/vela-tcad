#!/usr/bin/env python3
"""Replay one TransportModels failure with reference/precision A/B diagnostics.

The script first regenerates the failed Newton state with the extended restart
format, then produces a nearest-contact electron quasi-Fermi reference field.
Both coordinate representations are replayed in frozen-state mode.  Terminal
CSV output contains naive double, Neumaier-compensated, and long-double SG
reference totals computed from exactly the same saved state.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, getcontext
import json
import math
from pathlib import Path
import subprocess
import sys


getcontext().prec = 50


def run_runner(runner: Path, config: Path, allow_failure: bool = False) -> int:
    completed = subprocess.run([str(runner), "--config", str(config)], check=False)
    if completed.returncode != 0 and not allow_failure:
        raise RuntimeError(
            f"runner failed with exit code {completed.returncode}: {config}")
    return completed.returncode


def redirect_diagnostics(config: dict, output_dir: Path) -> None:
    sweep = config["sweep"]
    sweep["write_state_file"] = str(output_dir / "final_state.csv")
    sweep["write_state_every_point_prefix"] = str(output_dir / "state")
    diagnostics = sweep.setdefault("diagnostics", {})
    for name, filename in (
        ("terminal_balance", "terminal_balance.csv"),
        ("srh_balance", "srh_balance.csv"),
        ("contact_edge", "contact_edges.csv"),
        ("newton_history", "newton_history.csv"),
    ):
        section = diagnostics.setdefault(name, {})
        section["enabled"] = True
        section["csv_file"] = str(output_dir / filename)
    diagnostics["newton_history"]["rejected_state_directory"] = str(
        output_dir / "rejected_states")
    config["output_csv"] = str(output_dir / "curve.csv")


def require_extended_state(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as stream:
        header = next(csv.reader(stream))
    required = {
        "electron_qf_increment_V",
        "hole_qf_increment_V",
        "electron_qf_reference_V",
        "hole_qf_reference_V",
    }
    missing = sorted(required.difference(header))
    if missing:
        raise RuntimeError(
            f"state does not contain cancellation-free coordinates: {missing}")
    return header


def contact_map(mesh: dict) -> dict[str, list[int]]:
    return {entry["name"]: entry["node_ids"] for entry in mesh["contacts"]}


def build_partitioned_state(
    source: Path,
    destination: Path,
    mesh_path: Path,
    electron_contact_biases: dict[str, float],
) -> dict:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    contacts = contact_map(mesh)
    seeds: list[tuple[int, str, float]] = []
    for name, bias in electron_contact_biases.items():
        if name not in contacts:
            raise RuntimeError(f"mesh has no contact named {name!r}")
        seeds.extend((int(node), name, bias) for node in contacts[name])

    basin: dict[int, tuple[str, float]] = {}
    counts = {name: 0 for name in electron_contact_biases}
    for node_id, (x, y) in coordinates.items():
        seed_node, name, bias = min(
            seeds,
            key=lambda item: (
                (x - coordinates[item[0]][0]) ** 2
                + (y - coordinates[item[0]][1]) ** 2,
                item[1],
                item[0],
            ),
        )
        del seed_node
        basin[node_id] = (name, bias)
        counts[name] += 1

    with source.open(newline="", encoding="utf-8") as input_stream:
        reader = csv.DictReader(input_stream)
        if reader.fieldnames is None:
            raise RuntimeError(f"empty restart state: {source}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    max_reconstruction_error = Decimal(0)
    for row in rows:
        node_id = int(row["node_id"])
        _, new_reference_float = basin[node_id]
        old_reference = Decimal(row["electron_qf_reference_V"])
        old_increment = Decimal(row["electron_qf_increment_V"])
        new_reference_text = format(new_reference_float, ".17g")
        new_reference = Decimal(new_reference_text)
        new_increment = old_reference + old_increment - new_reference
        row["electron_qf_reference_V"] = new_reference_text
        row["electron_qf_increment_V"] = format(float(new_increment), ".17g")
        reconstructed = new_reference + Decimal(row["electron_qf_increment_V"])
        physical = Decimal(row["phin"])
        max_reconstruction_error = max(
            max_reconstruction_error, abs(reconstructed - physical))

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as output_stream:
        writer = csv.DictWriter(output_stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return {
        "basin_node_counts": counts,
        "max_physical_qf_reconstruction_error_V": float(max_reconstruction_error),
    }


def make_frozen_config(base: dict, state: Path, output_dir: Path) -> dict:
    config = json.loads(json.dumps(base))
    config["_comment"] = "Frozen replay for quasi-Fermi reference/precision A/B"
    redirect_diagnostics(config, output_dir)
    config["solver"]["method"] = "frozen_state"
    sweep = config["sweep"]
    sweep["initial_state_file"] = str(state)
    sweep["frozen_state_compute_current"] = True
    sweep["write_state_file"] = str(output_dir / "replayed_state.csv")
    sweep.pop("write_state_every_point_prefix", None)
    sweep["max_retries"] = 0
    return config


def read_terminal(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {row["contact"]: row for row in csv.DictReader(stream)}


def finite_float(row: dict[str, str], column: str) -> float:
    value = float(row[column])
    if not math.isfinite(value):
        raise RuntimeError(f"non-finite {column}: {row[column]}")
    return value


def summarize(global_dir: Path, partition_dir: Path, output_dir: Path) -> list[dict]:
    global_rows = read_terminal(global_dir / "terminal_balance.csv")
    partition_rows = read_terminal(partition_dir / "terminal_balance.csv")
    columns = (
        "current_electron_A_per_um",
        "current_electron_compensated_A_per_um",
        "current_electron_long_double_reference_A_per_um",
        "current_hole_A_per_um",
        "current_total_A_per_um",
    )
    records: list[dict] = []
    for contact in sorted(global_rows.keys() & partition_rows.keys()):
        record: dict[str, object] = {"contact": contact}
        for label, rows in (("global", global_rows), ("partition", partition_rows)):
            for column in columns:
                record[f"{label}_{column}"] = finite_float(rows[contact], column)
        global_e = record["global_current_electron_A_per_um"]
        partition_e = record["partition_current_electron_A_per_um"]
        record["partition_minus_global_electron_A_per_um"] = partition_e - global_e
        record["global_compensated_minus_naive_electron_A_per_um"] = (
            record["global_current_electron_compensated_A_per_um"] - global_e
        )
        record["global_long_double_minus_naive_electron_A_per_um"] = (
            record["global_current_electron_long_double_reference_A_per_um"] - global_e
        )
        record["partition_compensated_minus_naive_electron_A_per_um"] = (
            record["partition_current_electron_compensated_A_per_um"] - partition_e
        )
        record["partition_long_double_minus_naive_electron_A_per_um"] = (
            record["partition_current_electron_long_double_reference_A_per_um"]
            - partition_e
        )
        records.append(record)

    summary_csv = output_dir / "summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reuse-global-state", action="store_true")
    args = parser.parse_args()

    runner = args.runner.resolve()
    base_path = args.base_config.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    base = json.loads(base_path.read_text(encoding="utf-8"))

    global_solve_dir = output_dir / "global_solve"
    global_solve_dir.mkdir(parents=True, exist_ok=True)
    global_config = json.loads(json.dumps(base))
    redirect_diagnostics(global_config, global_solve_dir)
    global_config_path = global_solve_dir / "config.json"
    global_config_path.write_text(
        json.dumps(global_config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    global_state = global_solve_dir / "final_state.csv"
    solve_return_code = 0
    if not args.reuse_global_state or not global_state.exists():
        solve_return_code = run_runner(runner, global_config_path, allow_failure=True)
    require_extended_state(global_state)

    contact_biases = {
        entry["name"]: float(entry["bias"])
        for entry in base["contacts"]
    }
    # TransportModels uses an n+ polysilicon gate in addition to n+ source
    # and drain.  Include every such electron-majority terminal so the
    # partition does not create an artificial reference jump at the gate.
    electron_contacts = {
        name: contact_biases[name]
        for name in ("source", "drain", "gate")
        if name in contact_biases
    }
    partition_state = output_dir / "partition_state.csv"
    partition_metadata = build_partitioned_state(
        global_state,
        partition_state,
        Path(base["mesh_file"]),
        electron_contacts,
    )

    for label, state in (("global_replay", global_state),
                         ("partition_replay", partition_state)):
        replay_dir = output_dir / label
        replay_dir.mkdir(parents=True, exist_ok=True)
        config = make_frozen_config(base, state, replay_dir)
        config_path = replay_dir / "config.json"
        config_path.write_text(
            json.dumps(config, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        run_runner(runner, config_path)

    records = summarize(
        output_dir / "global_replay",
        output_dir / "partition_replay",
        output_dir,
    )
    execution = {
        "runner": str(runner),
        "base_config": str(base_path),
        "global_solve_return_code": solve_return_code,
        "global_state": str(global_state),
        "partition_state": str(partition_state),
        "partition": partition_metadata,
        "records": records,
    }
    (output_dir / "execution.json").write_text(
        json.dumps(execution, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - command-line failure path
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
