#!/usr/bin/env python3
"""Run the TransportModels Vg=1 V, Vd=2 V Sentaurus Frozen-Q oracle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
DEFAULT_OUTPUT = BASELINE / "frozen_q_oracle_vg1_vd2_run01"
DEFAULT_REPORT_JSON = (
    REPO_ROOT / "docs/validation/transportmodels_dg_frozen_q_oracle_2026-08-20.json"
)
DEFAULT_REPORT_MD = (
    REPO_ROOT / "docs/validation/transportmodels_dg_frozen_q_oracle_2026-08-20.md"
)
DEFAULT_RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
BASE_CONFIG = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "00_dg_idvd_curve.json"
)
SENTAURUS_EXPORT = BASELINE / "generated/sim_fields/dg_idvd"
REFERENCE_CURVE = (
    BASELINE
    / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvd_reference.csv"
)
SELF_CONSISTENT_CURVE = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "dg_idvd_curve_comparison_candidate.csv"
)
SELF_CONSISTENT_STATE = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "dg_idvd_curve_state_bias_2p000000.csv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest.upper()


def current_at(path: Path, column: str, bias_V: float) -> float:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if abs(float(row["bias_V"]) - bias_V) <= 1.0e-12:
                return abs(float(row[column]))
    raise RuntimeError(f"Missing bias {bias_V:g} V in {path}")


def merge_sentaurus_q_into_vela_state(
    sentaurus_restart: Path, output: Path
) -> None:
    with sentaurus_restart.open("r", encoding="utf-8-sig", newline="") as stream:
        sentaurus_rows = {
            int(row["node_id"]): row for row in csv.DictReader(stream)
        }
    with SELF_CONSISTENT_STATE.open(
        "r", encoding="utf-8-sig", newline=""
    ) as stream:
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if "electron_quantum_potential_like_V" not in fieldnames:
        fieldnames.append("electron_quantum_potential_like_V")
    if len(rows) != len(sentaurus_rows):
        raise RuntimeError("Vela and Sentaurus restart node counts differ")
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            node_id = int(row["node_id"])
            row["electron_quantum_potential_V"] = sentaurus_rows[node_id][
                "electron_quantum_potential_V"
            ]
            row["electron_quantum_potential_like_V"] = sentaurus_rows[node_id][
                "electron_quantum_potential_like_V"
            ]
            writer.writerow(row)


def make_config(output_dir: Path, restart: Path) -> Path:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    config["_comment"] = (
        "Phase-1 TransportModels diagnostic: solve DD variables while holding "
        "the imported Sentaurus electron quantum potential fixed."
    )
    quantum = config["solver"]["electron_quantum_potential"]
    quantum["enabled"] = True
    quantum["coupling_mode"] = "frozen"
    config["solver"]["verbose"] = False
    for contact in config["contacts"]:
        if contact["name"].lower() == "drain":
            contact["bias"] = 2.0
    config["output_csv"] = str((output_dir / "frozen_q_idvd_2V.csv").resolve())
    config["log_file"] = str((output_dir / "frozen_q_idvd_2V.log").resolve())
    sweep = config["sweep"]
    sweep.update(
        {
            "start": 2.0,
            "stop": 2.0,
            "step": 0.1,
            "bias_points": [2.0],
            "initial_state_file": str(restart.resolve()),
            "write_vtk": False,
            "write_state_file": str((output_dir / "frozen_q_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((output_dir / "frozen_q_state").resolve()),
        }
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["transport"] = {"enabled": True}
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((output_dir / "frozen_q_terminal_balance.csv").resolve()),
    }
    config_path = output_dir / "frozen_q_idvd_2V.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def render_markdown(summary: dict[str, object]) -> str:
    current = summary["current_A_per_um"]
    error = summary["relative_error"]
    reduction = summary["error_reduction"]
    assert isinstance(current, dict) and isinstance(error, dict) and isinstance(reduction, dict)
    classification = summary["classification"]
    return f"""# TransportModels DG Frozen-Q oracle

Work point: Vg = 1.0 V, Vd = 2.0 V

Status: **{summary['status']}**

## Terminal current

| Result | Id (A/um) | Relative error versus Sentaurus |
|---|---:|---:|
| Sentaurus DG | {current['sentaurus']:.12g} | 0 |
| Vela self-consistent DG baseline | {current['vela_self_consistent']:.12g} | {error['vela_self_consistent']:.6%} |
| Vela with Sentaurus Frozen-Q | {current['vela_frozen_q']:.12g} | {error['vela_frozen_q']:.6%} |

Frozen-Q removes **{reduction['absolute_percentage_points']:.6f} percentage
points** or **{reduction['fraction_of_baseline_error_removed']:.2%}** of the
self-consistent endpoint current error.

## Interpretation

Classification: **{classification}**.

The Frozen-Q run changes only the DD variables while preserving the imported
Sentaurus electron quantum potential. The result therefore separates the DG
field/equation contribution from the classical transport and mobility path.
It is a diagnostic oracle, not a production configuration.

The converged Vela self-consistent 2 V state supplies the initial electrostatic
and carrier variables; only its electron quantum potential is replaced by the
Sentaurus value. This avoids attributing an initial-state representation
mismatch to the Frozen-Q experiment.

## Provenance

- Imported restart SHA-256: `{summary['hashes']['sentaurus_restart']}`
- Hybrid restart SHA-256: `{summary['hashes']['hybrid_restart']}`
- Config SHA-256: `{summary['hashes']['config']}`
- Final state SHA-256: `{summary['hashes']['final_state']}`
- Curve SHA-256: `{summary['hashes']['curve']}`
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT_JSON)
    parser.add_argument("--report-md", type=Path, default=DEFAULT_REPORT_MD)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    restart = output_dir / "sentaurus_dg_idvd_restart.csv"
    converter = REPO_ROOT / "scripts/sentaurus_fields_to_restart.py"
    subprocess.run(
        [
            sys.executable,
            str(converter),
            "--export-dir",
            str(SENTAURUS_EXPORT),
            "--output",
            str(restart),
        ],
        cwd=REPO_ROOT,
        check=True,
    )
    hybrid_restart = output_dir / "vela_state_with_sentaurus_q.csv"
    merge_sentaurus_q_into_vela_state(restart, hybrid_restart)
    config_path = make_config(output_dir, hybrid_restart)
    if not args.execute:
        print(config_path)
        return 0

    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(runner)
    console_path = output_dir / "frozen_q.console.log"
    process = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
    )
    console_path.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    if process.returncode != 0:
        print(console_path)
        return process.returncode

    curve_path = output_dir / "frozen_q_idvd_2V.csv"
    state_path = output_dir / "frozen_q_final_state.csv"
    sentaurus = current_at(REFERENCE_CURVE, "current_total", 2.0)
    self_consistent = current_at(SELF_CONSISTENT_CURVE, "current_total_A_per_um", 2.0)
    frozen_q = current_at(curve_path, "current_total_A_per_um", 2.0)
    self_error = abs(self_consistent - sentaurus) / sentaurus
    frozen_error = abs(frozen_q - sentaurus) / sentaurus
    removed = self_error - frozen_error
    if frozen_error <= 0.03:
        classification = "DG equation/field is the dominant endpoint error source"
    elif frozen_error < self_error:
        classification = "DG field contributes, but transport/mobility coupling remains material"
    else:
        classification = "Frozen-Q does not improve the endpoint; inspect transport/mobility first"

    summary: dict[str, object] = {
        "schema": "vela.transportmodels.dg_frozen_q_oracle.v1",
        "status": "pass",
        "work_point": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        "initialization": {
            "dd_variables": "converged Vela self-consistent DG state at Vg=1 V, Vd=2 V",
            "electron_quantum_potential": "Sentaurus DG state at Vg=1 V, Vd=2 V",
        },
        "current_A_per_um": {
            "sentaurus": sentaurus,
            "vela_self_consistent": self_consistent,
            "vela_frozen_q": frozen_q,
        },
        "relative_error": {
            "vela_self_consistent": self_error,
            "vela_frozen_q": frozen_error,
        },
        "error_reduction": {
            "absolute_percentage_points": 100.0 * removed,
            "fraction_of_baseline_error_removed": removed / self_error,
        },
        "classification": classification,
        "paths": {
            "sentaurus_restart": str(restart),
            "hybrid_restart": str(hybrid_restart),
            "config": str(config_path),
            "curve": str(curve_path),
            "final_state": str(state_path),
            "console": str(console_path),
        },
        "hashes": {
            "sentaurus_restart": sha256(restart),
            "hybrid_restart": sha256(hybrid_restart),
            "config": sha256(config_path),
            "curve": sha256(curve_path),
            "final_state": sha256(state_path),
        },
    }
    summary_path = output_dir / "frozen_q_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    markdown = render_markdown(summary)
    (output_dir / "frozen_q_summary.md").write_text(markdown, encoding="utf-8")
    report_json = args.report_json.resolve()
    report_md = args.report_md.resolve()
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    report_md.write_text(markdown, encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
