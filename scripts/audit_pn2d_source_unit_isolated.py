#!/usr/bin/env python3
"""Isolate the continuity source-integral factor on identical states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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


def run_config(
    runner: Path,
    source_config: Path,
    output_csv: Path,
    config_path: Path,
) -> list[dict[str, str]]:
    cfg = json.loads(source_config.read_text(encoding="utf-8"))
    cfg["output_csv"] = str(output_csv.resolve())
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(cfg, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if not output_csv.is_file():
        raise RuntimeError(
            f"{config_path.name} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    return rows(output_csv)


def relative(left: float, right: float) -> float:
    return abs(left - right) / max(1.0, abs(left), abs(right))


def ratio(after: float, before: float) -> float | None:
    if before == 0.0:
        return None
    return after / before


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor-one-runner", type=Path, required=True)
    parser.add_argument("--factor-scaled-runner", type=Path, required=True)
    parser.add_argument("--phase-e-root", type=Path, required=True)
    parser.add_argument("--task6-root", type=Path, required=True)
    parser.add_argument("--forward-iv-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    factor_one = args.factor_one_runner.resolve()
    factor_scaled = args.factor_scaled_runner.resolve()
    phase_e = args.phase_e_root.resolve()
    task6 = args.task6_root.resolve()
    forward_iv = args.forward_iv_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw"

    term_compare: list[dict[str, object]] = []
    edge_compare: list[dict[str, object]] = []
    jacobian_compare: list[dict[str, object]] = []
    update_compare: list[dict[str, object]] = []
    maximum_edge_difference = 0.0
    maximum_flux_term_difference = 0.0
    srh_ratios: list[float] = []
    impact_ratios: list[float] = []
    for topology in ("mirror", "sketch"):
        for bias in (1, 10, 20):
            tag = f"m{bias}V"
            source = phase_e / "raw" / topology / tag
            task6_source = task6 / "raw" / topology / tag
            work = raw / topology / tag
            probe_sources = {
                "terms": source / "vela_production_terms.json",
                "edges": source / "vela_production_edges.json",
                "jacobian": source / "jacobian.json",
                "carrier": source / "first_update.json",
                "coupled": (
                    task6_source
                    / "vela_global_qfp_config_coupled.json"
                ),
            }
            values: dict[
                tuple[str, str], list[dict[str, str]]
            ] = {}
            for runner_name, runner in (
                ("factor_one", factor_one),
                ("factor_scaled", factor_scaled),
            ):
                for probe, source_config in probe_sources.items():
                    values[(runner_name, probe)] = run_config(
                        runner,
                        source_config,
                        work / f"{runner_name}_{probe}.csv",
                        work / f"{runner_name}_{probe}.json",
                    )

            before_terms = values[("factor_one", "terms")]
            after_terms = values[("factor_scaled", "terms")]
            for node, (before, after) in enumerate(
                zip(before_terms, after_terms)
            ):
                for carrier in ("electron", "hole"):
                    flux_before = float(before[f"{carrier}_flux"])
                    flux_after = float(after[f"{carrier}_flux"])
                    srh_before = float(
                        before[f"{carrier}_recombination"]
                    )
                    srh_after = float(
                        after[f"{carrier}_recombination"]
                    )
                    impact_before = float(before[f"{carrier}_impact"])
                    impact_after = float(after[f"{carrier}_impact"])
                    srh_ratio = ratio(srh_after, srh_before)
                    impact_ratio = ratio(impact_after, impact_before)
                    if srh_ratio is not None:
                        srh_ratios.append(srh_ratio)
                    if impact_ratio is not None:
                        impact_ratios.append(impact_ratio)
                    maximum_flux_term_difference = max(
                        maximum_flux_term_difference,
                        relative(flux_before, flux_after),
                    )
                    term_compare.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "node_id": node,
                            "carrier": carrier,
                            "factor_one_flux": flux_before,
                            "factor_scaled_flux": flux_after,
                            "flux_relative_difference": relative(
                                flux_before, flux_after
                            ),
                            "factor_one_srh": srh_before,
                            "factor_scaled_srh": srh_after,
                            "srh_scaled_over_one": srh_ratio,
                            "factor_one_impact": impact_before,
                            "factor_scaled_impact": impact_after,
                            "impact_scaled_over_one": impact_ratio,
                            "factor_one_residual": float(
                                before[f"{carrier}_residual"]
                            ),
                            "factor_scaled_residual": float(
                                after[f"{carrier}_residual"]
                            ),
                        }
                    )

            before_edges = values[("factor_one", "edges")]
            after_edges = values[("factor_scaled", "edges")]
            for before, after in zip(before_edges, after_edges):
                for carrier in ("electron", "hole"):
                    flux_before = float(before[f"{carrier}_flux"])
                    flux_after = float(after[f"{carrier}_flux"])
                    difference = relative(flux_before, flux_after)
                    maximum_edge_difference = max(
                        maximum_edge_difference, difference
                    )
                    edge_compare.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "edge_id": before["edge_id"],
                            "carrier": carrier,
                            "factor_one_flux": flux_before,
                            "factor_scaled_flux": flux_after,
                            "relative_difference": difference,
                        }
                    )

            before_jac = {
                row["block"]: row
                for row in values[("factor_one", "jacobian")]
            }
            after_jac = {
                row["block"]: row
                for row in values[("factor_scaled", "jacobian")]
            }
            for block in sorted(before_jac):
                before = before_jac[block]
                after = after_jac[block]
                jacobian_compare.append(
                    {
                        "topology": topology,
                        "bias_V": -bias,
                        "block": block,
                        "factor_one_analytic_norm": before[
                            "analytic_norm"
                        ],
                        "factor_scaled_analytic_norm": after[
                            "analytic_norm"
                        ],
                        "scaled_over_one_analytic_norm": ratio(
                            float(after["analytic_norm"]),
                            float(before["analytic_norm"]),
                        ),
                        "factor_one_rel_diff": before["rel_diff"],
                        "factor_scaled_rel_diff": after["rel_diff"],
                    }
                )

            for mode in ("carrier", "coupled"):
                before_update = values[("factor_one", mode)]
                after_update = values[("factor_scaled", mode)]
                for before, after in zip(before_update, after_update):
                    node = int(before["node_id"])
                    if node not in (1, 5):
                        continue
                    for carrier, delta in (
                        ("electron", "delta_phin_V"),
                        ("hole", "delta_phip_V"),
                    ):
                        update_compare.append(
                            {
                                "topology": topology,
                                "bias_V": -bias,
                                "mode": (
                                    "carrier_only"
                                    if mode == "carrier"
                                    else "coupled"
                                ),
                                "node_id": node,
                                "carrier": carrier,
                                "factor_one_delta_qfp_V": before[delta],
                                "factor_scaled_delta_qfp_V": after[delta],
                                "scaled_over_one_abs_delta": (
                                    abs(float(after[delta]))
                                    / max(
                                        abs(float(before[delta])),
                                        1.0e-300,
                                    )
                                ),
                            }
                        )

    write_csv(output / "term_comparison.csv", term_compare)
    write_csv(output / "edge_flux_control.csv", edge_compare)
    write_csv(output / "jacobian_comparison.csv", jacobian_compare)
    write_csv(output / "first_update_comparison.csv", update_compare)
    finite_srh = [
        value for value in srh_ratios if math.isfinite(value)
    ]
    finite_impact = [
        value for value in impact_ratios if math.isfinite(value)
    ]
    maximum_srh_ratio_error = max(
        abs(value - 1.0e-8) for value in finite_srh
    )
    maximum_impact_ratio_error = (
        max(abs(value - 1.0e-8) for value in finite_impact)
        if finite_impact
        else 0.0
    )
    forward_manifest = json.loads(
        (forward_iv / "manifest.json").read_text(encoding="utf-8")
    )
    patch_rows = rows(forward_iv / "patch_effect.csv")
    maximum_forward_iv_relative_change = max(
        abs(float(row["post_over_pre_current"]) - 1.0)
        for row in patch_rows
        if row["deck"] == "unit_scaling"
        and row["post_over_pre_current"] != ""
    )
    dimensional_factor = (1.0e-6**2) / 1.0e-4
    outcome = (
        "source_factor_dimensionally_required_forward_iv_insensitive"
        if maximum_edge_difference <= 1.0e-14
        and maximum_flux_term_difference <= 1.0e-14
        and maximum_srh_ratio_error <= 1.0e-14
        and maximum_impact_ratio_error <= 1.0e-14
        and abs(dimensional_factor - 1.0e-8) <= 1.0e-20
        else "source_factor_audit_failed"
    )
    report = [
        "# Isolated continuity source-unit audit",
        "",
        f"Typed outcome: `{outcome}`.",
        "",
        "`factor_one` and `factor_scaled` are built from the same current "
        "source tree. The only changed expression is "
        "`continuitySourceIntegralFactor()`: `1` versus "
        "`(1e-6 m)^2/(1e-4 m^2/V/s) = 1e-8`.",
        "",
        f"Maximum SG edge-flux change: `{maximum_edge_difference:.6e}`.",
        "",
        f"Maximum carrier-term SG change: "
        f"`{maximum_flux_term_difference:.6e}`.",
        "",
        f"Maximum SRH ratio error from 1e-8: "
        f"`{maximum_srh_ratio_error:.6e}`.",
        "",
        f"Maximum impact ratio error from 1e-8: "
        f"`{maximum_impact_ratio_error:.6e}`.",
        "",
        f"The forward-IV smoke deck changes by at most "
        f"`{maximum_forward_iv_relative_change:.6e}` relative, so its prior "
        "agreement is not evidence that the source factor cancels.",
        "",
        "The subsequent division by the common continuity normalization "
        "`C0*D0` acts on both SG and source terms. It therefore does not cancel "
        "the relative area/mobility conversion needed before the terms are "
        "combined.",
        "",
    ]
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )
    outputs = (
        "term_comparison.csv",
        "edge_flux_control.csv",
        "jacobian_comparison.csv",
        "first_update_comparison.csv",
        "report.md",
    )
    manifest = {
        "schema_version": 1,
        "status": "valid" if "failed" not in outcome else "failed",
        "experiment": "pn2d_source_unit_isolated",
        "typed_outcome": outcome,
        "dimensional_source_integral_factor": dimensional_factor,
        "maximum_edge_flux_relative_difference": maximum_edge_difference,
        "maximum_flux_term_relative_difference": (
            maximum_flux_term_difference
        ),
        "maximum_srh_ratio_error_from_1e_8": maximum_srh_ratio_error,
        "maximum_impact_ratio_error_from_1e_8": (
            maximum_impact_ratio_error
        ),
        "maximum_forward_iv_relative_change": (
            maximum_forward_iv_relative_change
        ),
        "forward_iv_original_outcome": forward_manifest["typed_outcome"],
        "inputs": {
            "factor_one_runner": sha256(factor_one),
            "factor_scaled_runner": sha256(factor_scaled),
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
    return 0 if manifest["status"] == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())
