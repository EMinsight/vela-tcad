#!/usr/bin/env python3
"""Run fixed-Sentaurus-state A/B probes for the DG deep-off SRH gap."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
from pathlib import Path
from typing import Any, Callable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
FORMULA_SCRIPT = REPO / "scripts/run_transportmodels_sentaurus_formula_replay.py"
FORMULA_REPORT = REPO / "docs/validation/transportmodels_sentaurus_formula_replay_2026-08-23.json"
STRICT_ROOT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/scaled_filter"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/srh_ab"
REPORT = REPO / "docs/validation/transportmodels_dg_deep_off_srh_ab_2026-08-23.json"
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


def identity(_: dict[str, Any]) -> None:
    return None


def disable_fermi_bgn(config: dict[str, Any]) -> None:
    config["solver"]["bandgap_narrowing"]["fermi_statistics_correction"] = False


def use_boltzmann(config: dict[str, Any]) -> None:
    config["solver"]["carrier_statistics"]["model"] = "boltzmann"


def use_boltzmann_no_fermi_bgn(config: dict[str, Any]) -> None:
    use_boltzmann(config)
    disable_fermi_bgn(config)


def use_net_doping(config: dict[str, Any]) -> None:
    config["solver"]["srh_doping_dependence"]["concentration_basis"] = "net_doping"


def disable_doping_dependence(config: dict[str, Any]) -> None:
    config["solver"]["srh_doping_dependence"]["enabled"] = False


VARIANTS: tuple[tuple[str, Callable[[dict[str, Any]], None]], ...] = (
    ("baseline_fermi_oldslotboom", identity),
    ("no_fermi_bgn_correction", disable_fermi_bgn),
    ("boltzmann", use_boltzmann),
    ("boltzmann_no_fermi_bgn_correction", use_boltzmann_no_fermi_bgn),
    ("net_doping_lifetime", use_net_doping),
    ("constant_lifetime", disable_doping_dependence),
)


def run_variant(
    formula: Any,
    case: dict[str, Any],
    state: Path,
    name: str,
    mutate: Callable[[dict[str, Any]], None],
    run_dir: Path,
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    output_csv = run_dir / "carrier_terms.csv"
    config = formula.make_probe_config(
        case, "newton_carrier_term_probe", state, output_csv
    )
    mutate(config)
    config_path = run_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(run_dir / "probe.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=formula.runner_environment(),
        check=False,
    )
    (run_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"{name} at {case['bias_V']} V failed: {completed.stderr or completed.stdout}")
    return output_csv


def main() -> int:
    formula = load_formula_module()
    prior = json.loads(FORMULA_REPORT.read_text(encoding="utf-8"))
    prior_cases = {
        round(float(row["bias_V"]), 12): row
        for row in prior["cases"]
        if row["group"] == "dg_idvg_deep_off"
    }
    cases = {
        round(float(row["bias_V"]), 12): row
        for row in formula.cases()
        if row["group"] == "dg_idvg_deep_off"
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for bias in (-1.0, -0.84, -0.68):
        key = round(bias, 12)
        strict_balance = read_csv(STRICT_ROOT / tag(bias) / "srh_balance.csv")[-1]
        strict_generation = float(strict_balance["srh_generation_current_A_per_um"])
        strict_terms = read_csv(
            REF
            / "reports/transportmodels_dg_deep_off_strict_20260823/srh_decomposition"
            / tag(bias)
            / "newton_carrier_term.csv"
        )
        strict_sum = abs(sum(float(row["electron_recombination"]) for row in strict_terms))
        internal_to_current = strict_generation / strict_sum
        prior_case = prior_cases[key]
        sentaurus_srh = abs(
            float(prior_case["srh"]["sentaurus_signed_area_weighted_sum_cm-1_s-1"])
        ) * Q * 1.0e-12
        sentaurus_state = Path(prior_case["artifacts"]["sentaurus_state"])
        for variant, mutate in VARIANTS:
            term_csv = run_variant(
                formula,
                cases[key],
                sentaurus_state,
                variant,
                mutate,
                OUTPUT / tag(bias) / variant,
            )
            current = abs(
                sum(float(row["electron_recombination"]) for row in read_csv(term_csv))
            ) * internal_to_current
            rows.append(
                {
                    "bias_V": bias,
                    "variant": variant,
                    "sentaurus_exported_srh_A_per_um": sentaurus_srh,
                    "vela_formula_srh_A_per_um": current,
                    "relative_gap_fraction": abs(sentaurus_srh - current) / sentaurus_srh,
                    "signed_ratio_vela_over_sentaurus": current / sentaurus_srh,
                    "carrier_terms_csv": str(term_csv.resolve()),
                }
            )

    csv_path = OUTPUT / "srh_ab_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "vela.transportmodels.dg_deep_off_srh_ab.v1",
        "status": "complete",
        "method": (
            "Each A/B variant evaluates Vela production SRH terms on the same immutable "
            "Sentaurus electrostatic, quasi-Fermi, quantum-potential and density state. "
            "The mesh and barycentric source integration are unchanged."
        ),
        "rows": rows,
        "artifacts": {"csv": str(csv_path.resolve())},
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
