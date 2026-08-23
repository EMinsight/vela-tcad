#!/usr/bin/env python3
"""Seed the SLOT-LDMOS Stage 06 load line from the completed Stage 05 row."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


DEFAULT_STAGE05_IV = Path("outputs/stages/05_avalanche_on_60v/iv.csv")
DEFAULT_STAGE06_CONFIG = Path("simulation_06_bvds_external_resistor_final.json")
DEFAULT_OUTPUT_CONFIG = Path(
    "simulation_06_bvds_external_resistor_final_from_stage05.json"
)


class RestartPreparationError(ValueError):
    """Raised when Stage 05 cannot provide a consistent Stage 06 restart."""


def read_last_converged_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RestartPreparationError(f"Stage 05 IV file is empty: {path}")
    converged = [row for row in rows if row.get("converged") == "1"]
    if not converged:
        raise RestartPreparationError(
            f"Stage 05 IV file contains no converged row: {path}"
        )
    return converged[-1]


def prepare_restart_document(
    document: dict[str, Any],
    stage05_row: dict[str, str],
    max_inner_voltage_step_V: float,
    resume_boundary_control: bool = False,
) -> dict[str, Any]:
    if not math.isfinite(max_inner_voltage_step_V) or max_inner_voltage_step_V <= 0:
        raise RestartPreparationError(
            "max inner-voltage step must be finite and positive"
        )
    try:
        inner_voltage = float(stage05_row["inner_voltage_V"])
        outer_voltage = float(stage05_row["outer_voltage_V"])
    except (KeyError, TypeError, ValueError) as error:
        raise RestartPreparationError(
            "Stage 05 IV row must contain numeric inner_voltage_V and outer_voltage_V"
        ) from error
    if not math.isfinite(inner_voltage) or not math.isfinite(outer_voltage):
        raise RestartPreparationError("Stage 05 terminal voltages must be finite")

    sweep = document.get("sweep")
    if not isinstance(sweep, dict):
        raise RestartPreparationError("Stage 06 config is missing sweep")
    points = sweep.get("bias_points")
    if not isinstance(points, list) or not points:
        raise RestartPreparationError("Stage 06 config has no outer bias points")
    first_outer = float(points[0])
    tolerance = 1.0e-9 * max(1.0, abs(first_outer), abs(outer_voltage))
    if abs(first_outer - outer_voltage) > tolerance:
        raise RestartPreparationError(
            "Stage 05 final outer voltage does not match the first Stage 06 point: "
            f"{outer_voltage} V versus {first_outer} V"
        )

    circuit = sweep.get("external_circuit")
    if not isinstance(circuit, dict) or circuit.get("mode") != "series_resistor":
        raise RestartPreparationError(
            "Stage 06 must use external_circuit.mode=series_resistor"
        )
    circuit["initial_inner_voltage_V"] = inner_voltage
    circuit["max_inner_voltage_step_V"] = max_inner_voltage_step_V

    control = sweep.get("boundary_control")
    if not isinstance(control, dict):
        raise RestartPreparationError("Stage 06 config is missing boundary_control")
    control["resume"] = resume_boundary_control

    document["_stage06_restart"] = {
        "source": DEFAULT_STAGE05_IV.as_posix(),
        "stage05_outer_voltage_V": outer_voltage,
        "stage05_inner_voltage_V": inner_voltage,
        "max_inner_voltage_step_V": max_inner_voltage_step_V,
        "resume_boundary_control": resume_boundary_control,
        "reason": (
            "Preserve voltage/state consistency at the Stage 05 to Stage 06 "
            "handoff and give the monotone root bracket sufficient reach."
        ),
    }
    return document


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--stage05-iv", type=Path, default=DEFAULT_STAGE05_IV)
    parser.add_argument("--stage06-config", type=Path, default=DEFAULT_STAGE06_CONFIG)
    parser.add_argument("--output-config", type=Path, default=DEFAULT_OUTPUT_CONFIG)
    parser.add_argument("--max-inner-voltage-step", type=float, default=0.25)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse matching boundary-control checkpoints from an interrupted run.",
    )
    return parser.parse_args()


def resolve_under_bundle(bundle: Path, value: Path) -> Path:
    return value if value.is_absolute() else bundle / value


def main() -> int:
    args = parse_args()
    bundle = args.bundle.resolve()
    iv_path = resolve_under_bundle(bundle, args.stage05_iv)
    config_path = resolve_under_bundle(bundle, args.stage06_config)
    output_path = resolve_under_bundle(bundle, args.output_config)

    row = read_last_converged_row(iv_path)
    with config_path.open(encoding="utf-8") as handle:
        document = json.load(handle)
    prepared = prepare_restart_document(
        document, row, args.max_inner_voltage_step, args.resume
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(prepared, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output_config": str(output_path),
                "initial_inner_voltage_V": prepared["sweep"]["external_circuit"][
                    "initial_inner_voltage_V"
                ],
                "max_inner_voltage_step_V": prepared["sweep"]["external_circuit"][
                    "max_inner_voltage_step_V"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
