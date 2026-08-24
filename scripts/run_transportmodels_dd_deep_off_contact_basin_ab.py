#!/usr/bin/env python3
"""Run a three-point DD deep-off A/B with contact-basin QF references."""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
BASE_CONFIG = (
    REF
    / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24/runs/dd"
    / "03_dd_idvg_curve.json"
)
BASE_STATE = BASE_CONFIG.with_name("dd_idvg_final_bias_relax_final_state.csv")
BASE_CURVE = (
    REF
    / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24"
    / "dd_idvg_completed_prefix.csv"
)
SENT_CURVE = REF / "run02/normalized/dd_idvg.csv"
CONTACTS = (
    REF
    / "sentaurus_vm_runs/dd_deep_off_spatial_oracles_20260824/exports/vg_m1p00"
    / "contacts.csv"
)
OUTPUT = REF / "reports/transportmodels_dd_deep_off_contact_basin_ab_20260824"
RUNNER = REPO / "build-release/vela_example_runner.exe"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def nearest(rows: list[dict[str, str]], bias: float) -> dict[str, str]:
    return min(rows, key=lambda row: abs(float(row["bias_V"]) - bias))


def state_path(bias: float) -> Path:
    sign = "m" if bias < 0 else "p"
    tag = f"{sign}{abs(bias):.6f}".replace(".", "p")
    return OUTPUT / f"state_bias_{tag}.csv"


def contact_ids(name: str) -> set[int]:
    row = next(row for row in read_csv(CONTACTS) if row["name"] == name)
    return {int(value) for value in row["node_ids"].split(";") if value}


def reference_audit(path: Path, nodes: set[int]) -> dict[str, Any]:
    selected = [row for row in read_csv(path) if int(row["node_id"]) in nodes]
    electron = sorted({float(row["electron_qf_reference_V"]) for row in selected})
    hole = sorted({float(row["hole_qf_reference_V"]) for row in selected})
    return {
        "node_count": len(selected),
        "electron_reference_values_V": electron,
        "hole_reference_values_V": hole,
    }


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    config["_comment"] = (
        "DD deep-off causal A/B: fixed contract plus contact_basin quasi-Fermi references"
    )
    config["solver"]["quasi_fermi_reference"] = "contact_basin"
    config["output_csv"] = str((OUTPUT / "curve.csv").resolve())
    sweep = config["sweep"]
    sweep["start"] = -1.0
    sweep["stop"] = -0.68
    sweep["step"] = 0.16
    sweep["bias_points"] = [-1.0, -0.84, -0.68]
    sweep["initial_state_file"] = str(BASE_STATE.resolve())
    sweep["write_state_file"] = str((OUTPUT / "final_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str((OUTPUT / "state").resolve())
    diagnostics = sweep["diagnostics"]
    diagnostics["terminal_balance"]["csv_file"] = str(
        (OUTPUT / "terminal_balance.csv").resolve()
    )
    diagnostics["srh_balance"]["csv_file"] = str(
        (OUTPUT / "srh_balance.csv").resolve()
    )
    config_path = OUTPUT / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    env = os.environ.copy()
    env["Path"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("Path", "")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(OUTPUT / "runner.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )
    (OUTPUT / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (OUTPUT / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)

    sent = read_csv(SENT_CURVE)
    baseline = read_csv(BASE_CURVE)
    candidate = read_csv(OUTPUT / "curve.csv")
    balances = read_csv(OUTPUT / "srh_balance.csv")
    rows: list[dict[str, Any]] = []
    for bias in (-1.0, -0.84, -0.68):
        sent_row = nearest(sent, bias)
        base_row = nearest(baseline, bias)
        cand_row = nearest(candidate, bias)
        balance = nearest(balances, bias)
        sent_id = float(sent_row["current_total"])
        base_id = float(base_row["current_total_A_per_um"])
        cand_id = float(cand_row["current_total_A_per_um"])
        rows.append(
            {
                "gate_bias_V": bias,
                "sentaurus_Id_A_per_um": sent_id,
                "baseline_none_Id_A_per_um": base_id,
                "contact_basin_Id_A_per_um": cand_id,
                "baseline_relative_error": abs(base_id - sent_id) / abs(sent_id),
                "contact_basin_relative_error": abs(cand_id - sent_id) / abs(sent_id),
                "contact_basin_kcl_residual_A_per_um": float(
                    balance["four_terminal_kcl_residual_A_per_um"]
                ),
                "contact_basin_id_to_kcl_ratio": float(
                    balance["id_to_kcl_residual_ratio"]
                ),
                "contact_basin_srh_A_per_um": float(
                    balance["srh_net_current_A_per_um"]
                ),
            }
        )
    write_csv(OUTPUT / "comparison.csv", rows)
    report = {
        "schema": "vela.transportmodels_dd_deep_off_contact_basin_ab.v1",
        "status": "complete",
        "rows": rows,
        "reference_audit": {
            "source": reference_audit(state_path(-0.68), contact_ids("source")),
            "drain": reference_audit(state_path(-0.68), contact_ids("drain")),
        },
        "config": str(config_path.resolve()),
    }
    (OUTPUT / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "complete", "points": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
