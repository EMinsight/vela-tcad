#!/usr/bin/env python3
"""Compare PN2D avalanche-off BV mobility doping concentration bases."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


BASES = (
    "net_doping",
    "total_impurity",
    "cell_reconstructed_total_impurity",
)
COLORS = {
    "Sentaurus avalanche-off": "#5B5B5B",
    "net_doping": "#4C78A8",
    "total_impurity": "#F2A541",
    "cell_reconstructed_total_impurity": "#E45756",
}
MARKERS = {
    "Sentaurus avalanche-off": "s",
    "net_doping": "o",
    "total_impurity": "^",
    "cell_reconstructed_total_impurity": "D",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sentaurus_off_by_bias(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_csv(path):
        bias = int(round(abs(float(row["bias_V"]))))
        if "sentaurus_avalanche_off_A_per_um" in row:
            value = row["sentaurus_avalanche_off_A_per_um"]
        else:
            value = row["current_total_A_per_um"]
        result[bias] = abs(float(value))
    return result


def configure_case(
    base: dict[str, Any],
    basis: str,
    case_dir: Path,
) -> tuple[dict[str, Any], Path]:
    cfg = json.loads(json.dumps(base))
    cfg["solver"]["mobility"]["doping_concentration_basis"] = basis
    cfg["solver"]["impact_ionization"] = {"model": "none"}
    cfg["output_csv"] = str((case_dir / f"{basis}.csv").resolve())
    sweep = cfg["sweep"]
    sweep["start"] = 0.0
    sweep["stop"] = -20.0
    sweep["bias_points"] = [float(-value) for value in range(21)]
    sweep["write_vtk"] = False
    sweep["vtk_prefix"] = str((case_dir / basis).resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((case_dir / f"{basis}_newton.csv").resolve()),
    }
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((case_dir / f"{basis}_terminal_balance.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((case_dir / f"{basis}_continuity_balance.csv").resolve()),
    }
    diagnostics["terminal_current_method_compare"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((case_dir / f"{basis}_terminal_method.csv").resolve()),
    }
    config_path = case_dir / f"simulation_{basis}.json"
    return cfg, config_path


def run_case(
    runner: Path,
    base: dict[str, Any],
    basis: str,
    out_dir: Path,
) -> Path:
    case_dir = out_dir / basis
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg, config_path = configure_case(base, basis, case_dir)
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    (case_dir / "runner.stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (case_dir / "runner.stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"{basis} failed with exit code {completed.returncode}: "
            f"{completed.stderr[-2000:]}"
        )
    return Path(cfg["output_csv"])


def make_plot(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    series: dict[str, list[tuple[float, float]]] = {}
    sentaurus_seen: set[int] = set()
    for row in rows:
        basis = str(row["basis"])
        bias = abs(float(row["bias_V"]))
        series.setdefault(basis, []).append(
            (bias, float(row["vela_current_A_per_um"]))
        )
        ibias = int(round(bias))
        if ibias not in sentaurus_seen:
            series.setdefault("Sentaurus avalanche-off", []).append(
                (bias, float(row["sentaurus_current_A_per_um"]))
            )
            sentaurus_seen.add(ibias)

    fig, axis = plt.subplots(figsize=(9.2, 5.4), dpi=180)
    for name, points in series.items():
        points.sort()
        axis.semilogy(
            [point[0] for point in points],
            [max(point[1], 1.0e-30) for point in points],
            color=COLORS[name],
            marker=MARKERS[name],
            markersize=4.2,
            linewidth=2.0 if name == "Sentaurus avalanche-off" else 1.7,
            linestyle="--" if name == "Sentaurus avalanche-off" else "-",
            label=name,
        )
    axis.set_title("PN2D avalanche-off reverse current by mobility doping basis")
    axis.set_xlabel("Reverse bias |V| (V)")
    axis.set_ylabel("|Anode current| (A/um)")
    axis.set_xlim(0.0, 20.0)
    axis.grid(True, which="both", color="#D7D7D7", linewidth=0.6, alpha=0.75)
    axis.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-control", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sentaurus = sentaurus_off_by_bias(args.sentaurus_control)

    comparison: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "comparison": "PN2D coarse7x3 avalanche-off reverse BV",
        "biases_V": list(range(21)),
        "bases": {},
        "pairwise": {},
    }
    basis_currents: dict[str, dict[int, float]] = {}
    for basis in BASES:
        output = run_case(args.runner.resolve(), base, basis, args.out_dir.resolve())
        rows = read_csv(output)
        converged_rows = [row for row in rows if row["converged"] == "1"]
        current = {
            int(round(abs(float(row["bias_V"])))): abs(
                float(row["current_total_A_per_um"])
            )
            for row in converged_rows
        }
        basis_currents[basis] = current
        log_errors: list[float] = []
        relative_errors: list[float] = []
        closure_ratios: list[float] = []
        for bias in range(1, 21):
            vela = current[bias]
            reference = sentaurus[bias]
            log_error = math.log10(vela / reference)
            relative_error = (vela - reference) / reference
            log_errors.append(log_error)
            relative_errors.append(abs(relative_error))
            row = next(
                item
                for item in converged_rows
                if int(round(abs(float(item["bias_V"])))) == bias
            )
            closure = max(
                float(row["global_electron_continuity_closure_ratio"]),
                float(row["global_hole_continuity_closure_ratio"]),
            )
            closure_ratios.append(closure)
            comparison.append(
                {
                    "basis": basis,
                    "bias_V": -bias,
                    "vela_current_A_per_um": vela,
                    "sentaurus_current_A_per_um": reference,
                    "vela_over_sentaurus": vela / reference,
                    "signed_relative_error": relative_error,
                    "log10_ratio": log_error,
                    "global_closure_max_ratio": closure,
                    "converged": 1,
                }
            )
        summary["bases"][basis] = {
            "converged_points": len(converged_rows),
            "expected_points": 21,
            "log10_ratio_rmse": math.sqrt(
                sum(value * value for value in log_errors) / len(log_errors)
            ),
            "median_absolute_relative_error": sorted(relative_errors)[
                len(relative_errors) // 2
            ],
            "max_global_closure_ratio": max(closure_ratios),
            "ratio_at_minus_1V": current[1] / sentaurus[1],
            "ratio_at_minus_5V": current[5] / sentaurus[5],
            "ratio_at_minus_10V": current[10] / sentaurus[10],
            "ratio_at_minus_15V": current[15] / sentaurus[15],
            "ratio_at_minus_20V": current[20] / sentaurus[20],
        }

    baseline = basis_currents["net_doping"]
    for basis in BASES[1:]:
        deltas = {
            bias: (
                basis_currents[basis][bias] - baseline[bias]
            ) / baseline[bias]
            for bias in range(1, 21)
        }
        max_bias = max(deltas, key=lambda bias: abs(deltas[bias]))
        summary["pairwise"][f"{basis}_vs_net_doping"] = {
            "max_absolute_relative_difference": abs(deltas[max_bias]),
            "bias_at_max_difference_V": -max_bias,
            "signed_relative_difference_at_minus_1V": deltas[1],
            "signed_relative_difference_at_minus_10V": deltas[10],
            "signed_relative_difference_at_minus_20V": deltas[20],
        }

    write_csv(args.out_dir / "comparison.csv", comparison)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    make_plot(args.out_dir / "bv_off_doping_basis_compare.png", comparison)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
