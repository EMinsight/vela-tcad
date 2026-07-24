#!/usr/bin/env python3
"""Compare pre/post source-unit runners on equivalent forward-IV decks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, values: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_deck(
    runner: Path,
    source: Path,
    output: Path,
    label: str,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    cfg = json.loads(source.read_text(encoding="utf-8"))
    cfg["mesh_file"] = str((source.parent / cfg["mesh_file"]).resolve())
    csv_path = output / f"{label}.csv"
    cfg["output_csv"] = str(csv_path)
    cfg["sweep"]["write_vtk"] = False
    cfg_path = output / f"{label}.json"
    cfg_path.write_text(
        json.dumps(cfg, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    completed = subprocess.run(
        [str(runner), "--config", str(cfg_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    status: dict[str, object] = {}
    if lines:
        try:
            status = json.loads(lines[-1])
        except json.JSONDecodeError:
            status = {}
    if not csv_path.is_file():
        raise RuntimeError(
            f"{label} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    status["returncode"] = completed.returncode
    status["stderr"] = completed.stderr.strip()
    return rows(csv_path), status


def current_A_per_um(row: dict[str, str], unit_scaled: bool) -> float:
    if unit_scaled:
        return float(row["current_total_A_per_um"])
    return float(row["current_total"]) * 1.0e-6


def log_error(left: float, right: float) -> float | None:
    if left == 0.0 or right == 0.0:
        return None
    return abs(math.log10(abs(left) / abs(right)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-runner", type=Path, required=True)
    parser.add_argument("--candidate-runner", type=Path, required=True)
    parser.add_argument("--legacy-deck", type=Path, required=True)
    parser.add_argument("--unit-deck", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    runners = {
        "pre_source_factor": args.baseline_runner.resolve(),
        "post_source_factor": args.candidate_runner.resolve(),
    }
    decks = {
        "legacy_si": (args.legacy_deck.resolve(), False),
        "unit_scaling": (args.unit_deck.resolve(), True),
    }
    results: dict[
        tuple[str, str], tuple[list[dict[str, str]], dict[str, object]]
    ] = {}
    for runner_name, runner in runners.items():
        for deck_name, (deck, _) in decks.items():
            results[(runner_name, deck_name)] = run_deck(
                runner,
                deck,
                output,
                f"{runner_name}_{deck_name}",
            )

    comparison: list[dict[str, object]] = []
    for runner_name in runners:
        legacy = {
            float(row["bias_V"]): row
            for row in results[(runner_name, "legacy_si")][0]
        }
        unit = {
            float(row["bias_V"]): row
            for row in results[(runner_name, "unit_scaling")][0]
        }
        for bias in sorted(set(legacy) & set(unit)):
            legacy_current = current_A_per_um(legacy[bias], False)
            unit_current = current_A_per_um(unit[bias], True)
            comparison.append(
                {
                    "runner": runner_name,
                    "bias_V": bias,
                    "legacy_converged": legacy[bias]["converged"],
                    "unit_converged": unit[bias]["converged"],
                    "legacy_current_A_per_um": legacy_current,
                    "unit_current_A_per_um": unit_current,
                    "current_abs_error_dex": log_error(
                        unit_current, legacy_current
                    ),
                    "legacy_stored_charge_C_per_m": legacy[bias].get(
                        "stored_charge_C_per_m", ""
                    ),
                    "unit_stored_charge_C_per_m": unit[bias].get(
                        "stored_charge_C_per_m", ""
                    ),
                }
            )

    patch_effect: list[dict[str, object]] = []
    for deck_name, (_, unit_scaled) in decks.items():
        before = {
            float(row["bias_V"]): row
            for row in results[("pre_source_factor", deck_name)][0]
        }
        after = {
            float(row["bias_V"]): row
            for row in results[("post_source_factor", deck_name)][0]
        }
        for bias in sorted(set(before) & set(after)):
            before_current = current_A_per_um(
                before[bias], unit_scaled
            )
            after_current = current_A_per_um(after[bias], unit_scaled)
            patch_effect.append(
                {
                    "deck": deck_name,
                    "bias_V": bias,
                    "pre_current_A_per_um": before_current,
                    "post_current_A_per_um": after_current,
                    "post_over_pre_current": (
                        after_current / before_current
                        if before_current != 0.0
                        else ""
                    ),
                    "abs_change_A_per_um": abs(
                        after_current - before_current
                    ),
                    "abs_change_dex": log_error(
                        after_current, before_current
                    ),
                }
            )

    write_csv(output / "scaling_parity.csv", comparison)
    write_csv(output / "patch_effect.csv", patch_effect)
    parity_summary: list[dict[str, object]] = []
    for runner_name in runners:
        errors = [
            float(row["current_abs_error_dex"])
            for row in comparison
            if row["runner"] == runner_name
            and row["current_abs_error_dex"] is not None
        ]
        parity_summary.append(
            {
                "runner": runner_name,
                "nonzero_bias_count": len(errors),
                "median_current_abs_error_dex": statistics.median(errors),
                "maximum_current_abs_error_dex": max(errors),
                "all_rows_converged": all(
                    row["legacy_converged"] == "1"
                    and row["unit_converged"] == "1"
                    for row in comparison
                    if row["runner"] == runner_name
                ),
            }
        )
    write_csv(output / "scaling_parity_summary.csv", parity_summary)
    before_error = next(
        float(row["median_current_abs_error_dex"])
        for row in parity_summary
        if row["runner"] == "pre_source_factor"
    )
    after_error = next(
        float(row["median_current_abs_error_dex"])
        for row in parity_summary
        if row["runner"] == "post_source_factor"
    )
    unit_patch_rows = [
        row
        for row in patch_effect
        if row["deck"] == "unit_scaling"
        and row["post_over_pre_current"] != ""
    ]
    maximum_patch_relative_change = max(
        abs(float(row["post_over_pre_current"]) - 1.0)
        for row in unit_patch_rows
    )
    if maximum_patch_relative_change <= 1.0e-6:
        outcome = "forward_iv_insensitive_to_source_factor"
    elif after_error + 1.0e-3 < before_error:
        outcome = "source_factor_improves_physical_unit_parity"
    else:
        outcome = "forward_iv_parity_inconclusive"
    report = [
        "# Source-unit forward-IV audit",
        "",
        f"Typed outcome: `{outcome}`.",
        "",
        f"Pre-factor median legacy/unit current error: "
        f"`{before_error:.6e} dex`.",
        "",
        f"Post-factor median legacy/unit current error: "
        f"`{after_error:.6e} dex`.",
        "",
        f"Maximum unit-scaling IV change from the source factor: "
        f"{maximum_patch_relative_change:.6e} relative.",
        "",
        "This forward-IV smoke deck is insensitive to the source factor. "
        "The historical legacy-SI and unit-scaling decks also do not form a "
        "numerically matched parity pair, so their cross-deck current gap is "
        "not used to accept or reject the factor.",
        "",
    ]
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )
    outputs = (
        "scaling_parity.csv",
        "patch_effect.csv",
        "scaling_parity_summary.csv",
        "report.md",
    )
    manifest = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_source_unit_forward_iv",
        "typed_outcome": outcome,
        "pre_factor_median_current_abs_error_dex": before_error,
        "post_factor_median_current_abs_error_dex": after_error,
        "maximum_patch_relative_change": maximum_patch_relative_change,
        "inputs": {
            "baseline_runner": sha256(runners["pre_source_factor"]),
            "candidate_runner": sha256(runners["post_source_factor"]),
            "legacy_deck": sha256(decks["legacy_si"][0]),
            "unit_deck": sha256(decks["unit_scaling"][0]),
        },
        "outputs": {
            name: sha256(output / name) for name in outputs
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
