#!/usr/bin/env python3
"""Decompose the remaining hard-closed DG deep-off SRH current gap."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
FORMULA_SCRIPT = REPO / "scripts/run_transportmodels_sentaurus_formula_replay.py"
FORMULA_REPORT = REPO / "docs/validation/transportmodels_sentaurus_formula_replay_2026-08-23.json"
STRICT_ROOT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/scaled_filter"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/srh_decomposition"
REPORT = REPO / "docs/validation/transportmodels_dg_deep_off_strict_srh_2026-08-23.json"
RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19


def load_formula_module():
    spec = importlib.util.spec_from_file_location("transportmodels_formula_replay", FORMULA_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {FORMULA_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def tag(bias: float) -> str:
    return "m" + f"{abs(bias):.2f}".replace(".", "p")


def main() -> int:
    formula = load_formula_module()
    prior = json.loads(FORMULA_REPORT.read_text(encoding="utf-8"))
    prior_cases = {
        round(float(row["bias_V"]), 12): row
        for row in prior["cases"]
        if row["group"] == "dg_idvg_deep_off"
    }
    formula_cases = {
        round(float(row["bias_V"]), 12): row
        for row in formula.cases()
        if row["group"] == "dg_idvg_deep_off"
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for bias in (-1.0, -0.84, -0.68):
        key = round(bias, 12)
        run_dir = OUTPUT / tag(bias)
        run_dir.mkdir(parents=True, exist_ok=True)
        state = STRICT_ROOT / tag(bias) / "final_state.csv"
        term_csv, _ = formula.run_probe(
            RUNNER.resolve(),
            formula_cases[key],
            run_dir,
            "newton_carrier_term_probe",
            state,
        )
        strict_term_sum = sum(
            float(row["electron_recombination"]) for row in read_csv(term_csv)
        )
        balance = read_csv(STRICT_ROOT / tag(bias) / "srh_balance.csv")[-1]
        strict_generation = float(balance["srh_generation_current_A_per_um"])
        internal_to_A_per_um = strict_generation / abs(strict_term_sum)

        prior_case = prior_cases[key]
        sent_term_rows = read_csv(Path(prior_case["artifacts"]["carrier_terms"]))
        vela_formula_sent_state_sum = sum(
            float(row["electron_recombination"]) for row in sent_term_rows
        )
        vela_formula_sent_state_current = (
            abs(vela_formula_sent_state_sum) * internal_to_A_per_um
        )
        sent_export_integral = abs(
            float(prior_case["srh"]["sentaurus_signed_area_weighted_sum_cm-1_s-1"])
        ) * Q * 1.0e-12

        rows.append(
            {
                "bias_V": bias,
                "sentaurus_exported_srh_A_per_um": sent_export_integral,
                "vela_formula_on_sentaurus_state_A_per_um": vela_formula_sent_state_current,
                "vela_formula_on_strict_state_A_per_um": strict_generation,
                "sentaurus_to_vela_formula_gap_fraction": abs(
                    sent_export_integral - vela_formula_sent_state_current
                ) / sent_export_integral,
                "sentaurus_state_to_strict_state_gap_fraction": abs(
                    vela_formula_sent_state_current - strict_generation
                ) / vela_formula_sent_state_current,
                "full_sentaurus_to_strict_gap_fraction": abs(
                    sent_export_integral - strict_generation
                ) / sent_export_integral,
                "internal_residual_to_A_per_um": internal_to_A_per_um,
                "strict_carrier_terms": str(term_csv.resolve()),
                "sentaurus_carrier_terms": prior_case["artifacts"]["carrier_terms"],
            }
        )

    csv_path = OUTPUT / "strict_srh_decomposition.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "vela.transportmodels.dg_deep_off_strict_srh.v1",
        "status": "complete",
        "conversion_note": (
            "Sentaurus srhRecombination is integrated over um^2 mesh area and a 1 um "
            "device depth: q * sum(rate_cm-3 * area_um2) * 1e-12 gives A/um. "
            "Vela internal carrier-term sums are calibrated against the production "
            "strict-state SRH balance at the same bias."
        ),
        "points": rows,
        "artifacts": {"csv": str(csv_path.resolve())},
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
