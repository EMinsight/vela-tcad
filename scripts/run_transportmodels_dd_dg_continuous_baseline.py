#!/usr/bin/env python3
"""Run strict, single-pass TransportModels DD/DG baseline sweeps.

Each Id-Vg or Id-Vd curve is produced by one runner invocation with one
immutable configuration.  Adaptive internal voltage steps may be inserted,
but every trial starts only from the immediately preceding accepted state.
External restarts, historical same-bias seeds and pointwise reclosure are
intentionally unsupported.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import run_transportmodels_dd_dg_workflow as workflow
import transportmodels_fixed_contract as fixed


REPO = Path(__file__).resolve().parents[1]
DEFAULT_OVERLAY = (
    REPO / "configs/regression/transportmodels_dd_dg_continuous_v1.json"
)
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"


def load_overlay(path: Path = DEFAULT_OVERLAY) -> dict[str, Any]:
    payload = workflow.load_json(path)
    expected = "vela.transportmodels.dd_dg.continuous_numerical_contract.v1"
    if payload.get("schema") != expected:
        raise ValueError(f"Unsupported continuous-sweep contract: {path}")
    base = (path.parent / payload["base_contract"]).resolve()
    if fixed.load_contract(base)["contract_id"] != (
        "transportmodels-dd-dg-sentaurus2022-v1"
    ):
        raise ValueError(f"Unexpected base contract: {base}")
    return payload


def exact_curve_biases(branch: str, curve: str) -> list[float]:
    del branch  # DD and DG deliberately share the same bias lattice.
    contract = fixed.load_contract()
    key = "gate_bias_V" if curve == "idvg" else "drain_bias_V"
    return [float(value) for value in contract["bias_contract"][curve][key]]


def patch_manifest(
    manifest: dict[str, Any], run_dir: Path, overlay_path: Path = DEFAULT_OVERLAY
) -> dict[str, Any]:
    overlay = load_overlay(overlay_path)
    base_path = (overlay_path.parent / overlay["base_contract"]).resolve()

    def adaptive_sweep(sweep: dict[str, Any], start: float, stop: float,
                       nominal_step: float) -> None:
        sweep.pop("bias_points", None)
        sweep.update({
            "start": start,
            "stop": stop,
            "step": nominal_step,
            "initial_step": nominal_step,
            "min_step": overlay["adaptive_sweep"]["min_step_V"],
            "max_step": nominal_step,
            "shrink_factor": overlay["adaptive_sweep"]["shrink_factor"],
            "growth_factor": overlay["adaptive_sweep"]["growth_factor"],
            "max_retries": overlay["adaptive_sweep"]["max_retries"],
            "stop_on_failure": overlay["adaptive_sweep"]["stop_on_failure"],
        })

    for stage in manifest["stages"]:
        if "external_initial_state" in stage:
            raise ValueError(f"{stage['name']}: external restart is forbidden")
        path = Path(stage["config"])
        config = fixed.apply_contract(
            workflow.load_json(path), stage["branch"], base_path
        )
        fixed.deep_merge(config.setdefault("solver", {}), overlay["solver_numerics"])
        config["solver"]["verbose"] = False
        config["continuous_baseline_contract"] = {
            "id": overlay["contract_id"],
            "path": str(overlay_path.resolve()),
            "sha256": fixed.sha256(overlay_path),
            "rules": overlay["rules"],
        }
        if stage["name"].endswith(("_idvg_curve", "_idvd_curve")):
            curve = "idvg" if "_idvg_" in stage["name"] else "idvd"
            exact = exact_curve_biases(stage["branch"], curve)
            nominal_step = exact[1] - exact[0]
            sweep = config["sweep"]
            # The dependency stage already solved and saved exact[0].  Start
            # directly at exact[1], using that saved state as the immediately
            # preceding accepted point; re-solving exact[0] is both redundant
            # and can trigger line-search non-decrease at an identical bias.
            adaptive_sweep(sweep, exact[1], exact[-1], nominal_step)
            stage["execution_lattice"] = "adaptive_nominal_targets"
            stage["seed_bias_point"] = exact[0]
            stage["nominal_bias_points"] = exact[1:]
        elif stage["name"].endswith("_idvd_equilibrium"):
            init = overlay["idvd_initialization"][stage["branch"]]
            if init["mode"] == "continuous_gate_ramp":
                adaptive_sweep(
                    config["sweep"], init["gate_start_V"], init["gate_stop_V"],
                    init["gate_step_V"],
                )
                stage["execution_lattice"] = "adaptive_idvd_gate_initialization"
                stage["initialization_bias_points"] = [
                    init["gate_start_V"], init["gate_stop_V"]
                ]
            elif init["mode"] != "direct_equilibrium":
                raise ValueError(
                    f"{stage['name']}: unsupported initialization {init['mode']}"
                )
        violations = fixed.validate_config(config, stage["branch"], base_path)
        if violations:
            raise ValueError(f"{stage['name']}: fixed contract violations: {violations}")
        if config["solver"].get("quasi_fermi_reference") != "contact_basin":
            raise ValueError(f"{stage['name']}: contact_basin is not active")
        workflow.write_json(path, config)
        stage["config_sha256"] = workflow.sha256(path)

    audit_lineage(manifest, run_dir)
    manifest.update({
        "schema": "vela.transportmodels.dd_dg.continuous_baseline.v1",
        "continuous_contract": {
            "path": str(overlay_path.resolve()),
            "sha256": fixed.sha256(overlay_path),
            "rules": overlay["rules"],
        },
        "status": "materialized",
    })
    workflow.write_json(run_dir / "workflow_manifest.json", manifest)
    return manifest


def audit_lineage(manifest: dict[str, Any], run_dir: Path) -> None:
    stages = {stage["name"]: stage for stage in manifest["stages"]}
    root = run_dir.resolve()
    for stage in manifest["stages"]:
        config = workflow.load_json(Path(stage["config"]))
        initial = config["sweep"].get("initial_state_file")
        dependencies = stage.get("depends_on", [])
        if not dependencies:
            if initial is not None:
                raise ValueError(f"{stage['name']}: root stage has an initial state")
            continue
        if len(dependencies) != 1:
            raise ValueError(f"{stage['name']}: exactly one predecessor is required")
        predecessor = stages[dependencies[0]]
        if Path(initial).resolve() != Path(predecessor["final_state_file"]).resolve():
            raise ValueError(f"{stage['name']}: initial state is not its predecessor")
        if root not in Path(initial).resolve().parents:
            raise ValueError(f"{stage['name']}: initial state escapes the fresh run")
        if stage["name"].endswith(("_idvg_curve", "_idvd_curve")):
            if "bias_points" in config["sweep"]:
                raise ValueError(f"{stage['name']}: explicit bias_points disable adaptivity")


def write_summary(manifest: dict[str, Any], run_dir: Path) -> None:
    rows = []
    for stage in manifest["stages"]:
        rows.append({
            "name": stage["name"],
            "branch": stage["branch"],
            "status": stage.get("status", "not-run"),
            "returncode": stage.get("returncode"),
            "config_sha256": stage["config_sha256"],
            "predecessor": (stage.get("depends_on") or [None])[0],
        })
    summary = {
        "schema": "vela.transportmodels.dd_dg.continuous_baseline.summary.v1",
        "status": manifest.get("status"),
        "comparison_status": manifest.get("comparison_status"),
        "rules": manifest["continuous_contract"]["rules"],
        "stages": rows,
        "comparisons": manifest.get("comparisons", {}),
        "strict_acceptance": manifest.get("strict_acceptance"),
    }
    workflow.write_json(run_dir / "continuous_baseline_summary.json", summary)
    lines = [
        "# TransportModels DD/DG strict continuous baseline",
        "",
        f"- Solver status: `{summary['status']}`",
        f"- Comparison status: `{summary['comparison_status']}`",
        "- Curves use one immutable configuration and only the previous accepted state.",
        "- Adaptive internal steps are permitted; external restart, same-bias historical seed, pointwise reclosure and interpolation are forbidden.",
        "",
        "| Stage | Branch | Status | Predecessor | Config SHA-256 |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['branch']} | {row['status']} | "
            f"{row['predecessor'] or '-'} | `{row['config_sha256']}` |"
        )
    strict = summary.get("strict_acceptance") or {}
    if strict:
        lines.extend(["", "## Fixed-contract acceptance", ""])
        for branch, result in strict["branches"].items():
            idvg = result["idvg"]
            idvd = result["idvd"]
            lines.append(
                f"- {branch.upper()}: overall=`{result['pass']}`; "
                f"Id-Vg transition={idvg['transition_max_log_error_dex']:.6f} dex; "
                f"Id-Vg on={100 * idvg['on_max_relative_error']:.3f}%; "
                f"Id-Vd={100 * idvd['nonzero_max_relative_error']:.3f}%; "
                f"Id-Vd@2V={100 * idvd['endpoint_relative_error']:.3f}%"
            )
    (run_dir / "continuous_baseline_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def terminal_kcl_by_bias(
    manifest: dict[str, Any], curve_stage: dict[str, Any]
) -> dict[float, float]:
    stages = {stage["name"]: stage for stage in manifest["stages"]}
    seed_name = curve_stage["comparison_seed"]["stage"]
    result: dict[float, float] = {}
    for stage in (stages[seed_name], curve_stage):
        config = workflow.load_json(Path(stage["config"]))
        path = Path(
            config["sweep"]["diagnostics"]["terminal_balance"]["csv_file"]
        )
        grouped: dict[float, list[float]] = {}
        for row in read_csv(path):
            bias = round(float(row["bias_V"]), 12)
            grouped.setdefault(bias, []).append(
                float(row["current_total_A_per_um"])
            )
        result.update({bias: abs(math.fsum(values)) for bias, values in grouped.items()})
    return result


def aligned_curve(
    stage: dict[str, Any], manifest: dict[str, Any]
) -> list[dict[str, float]]:
    comparison = manifest["comparisons"][stage["name"]]
    reference_rows = read_csv(Path(stage["reference"]))
    candidate_rows = read_csv(Path(comparison["candidate"]))
    reference = {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in reference_rows
    }
    candidate = {
        round(float(row["bias_V"]), 12): abs(
            float(row["current_total_A_per_um"])
        )
        for row in candidate_rows
    }
    expected = [round(float(value), 12) for value in (
        stage["comparison_seed"]["reference_biases"]
    )]
    if sorted(reference) != expected or sorted(candidate) != expected:
        raise ValueError(f"{stage['name']}: reference/candidate lattice mismatch")
    return [{
        "bias_V": bias,
        "reference_A_per_um": reference[bias],
        "candidate_A_per_um": candidate[bias],
        "relative_error": abs(candidate[bias] - reference[bias])
        / max(reference[bias], 1.0e-300),
        "log_error_dex": abs(
            math.log10(max(candidate[bias], 1.0e-300))
            - math.log10(max(reference[bias], 1.0e-300))
        ),
    } for bias in expected]


def apply_strict_acceptance(
    manifest: dict[str, Any], run_dir: Path
) -> dict[str, Any]:
    limits = fixed.load_contract()["acceptance"]
    stages = {stage["name"]: stage for stage in manifest["stages"]}
    branches: dict[str, Any] = {}
    for branch in ("dd", "dg"):
        idvg_stage = stages[f"{branch}_idvg_curve"]
        idvd_stage = stages[f"{branch}_idvd_curve"]
        idvg = aligned_curve(idvg_stage, manifest)
        idvd = aligned_curve(idvd_stage, manifest)
        kcl = terminal_kcl_by_bias(manifest, idvg_stage)
        deep = []
        for row in idvg[:3]:
            residual = kcl.get(row["bias_V"], math.inf)
            ratio = row["candidate_A_per_um"] / residual if residual > 0 else math.inf
            deep.append({
                **row,
                "four_terminal_kcl_residual_A_per_um": residual,
                "id_to_kcl_ratio": ratio,
                "pass": (
                    row["log_error_dex"]
                    <= limits["idvg_transition_max_absolute_log_error_dex"]
                    and ratio >= limits["deep_off_min_id_to_kcl_ratio"]
                ),
            })
        transition = idvg[3:8]
        on = idvg[8:]
        nonzero_idvd = [row for row in idvd if row["bias_V"] > 0.0]
        idvg_metrics = {
            "deep_off": deep,
            "deep_off_pass": all(row["pass"] for row in deep),
            "transition_max_log_error_dex": max(
                row["log_error_dex"] for row in transition
            ),
            "on_max_relative_error": max(row["relative_error"] for row in on),
        }
        idvd_metrics = {
            "excluded_zero_bias_from_relative_error": True,
            "nonzero_points": len(nonzero_idvd),
            "nonzero_max_relative_error": max(
                row["relative_error"] for row in nonzero_idvd
            ),
            "endpoint_relative_error": idvd[-1]["relative_error"],
        }
        gates = {
            "deep_off": idvg_metrics["deep_off_pass"],
            "idvg_transition": idvg_metrics["transition_max_log_error_dex"]
            <= limits["idvg_transition_max_absolute_log_error_dex"],
            "idvg_on": idvg_metrics["on_max_relative_error"]
            <= limits["idvg_on_max_absolute_relative_error"],
            "idvd_nonzero": idvd_metrics["nonzero_max_relative_error"]
            <= limits["idvd_max_absolute_relative_error"],
            "idvd_endpoint": idvd_metrics["endpoint_relative_error"]
            <= limits["idvd_2V_max_absolute_relative_error"],
        }
        branches[branch] = {
            "pass": all(gates.values()),
            "gates": gates,
            "idvg": idvg_metrics,
            "idvd": idvd_metrics,
        }
    report = {
        "schema": "vela.transportmodels.dd_dg.strict_acceptance.v1",
        "limits": limits,
        "branches": branches,
        "overall_pass": all(result["pass"] for result in branches.values()),
    }
    workflow.write_json(run_dir / "strict_acceptance.json", report)
    manifest["strict_acceptance"] = report
    manifest["comparison_status"] = "pass" if report["overall_pass"] else "fail"
    workflow.write_json(run_dir / "workflow_manifest.json", manifest)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--overlay", type=Path, default=DEFAULT_OVERLAY)
    parser.add_argument("--postprocess-existing", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if args.postprocess_existing:
        manifest_path = output / "workflow_manifest.json"
        if not manifest_path.is_file():
            parser.error("--postprocess-existing requires workflow_manifest.json")
        manifest = workflow.load_json(manifest_path)
        expected_comparisons = {
            f"{branch}_{curve}_curve"
            for branch in ("dd", "dg") for curve in ("idvg", "idvd")
        }
        if set(manifest.get("comparisons", {})) != expected_comparisons:
            parser.error("existing manifest does not contain all four comparisons")
        apply_strict_acceptance(manifest, output)
        write_summary(manifest, output)
        print(json.dumps({
            "status": manifest.get("status"),
            "comparison_status": manifest.get("comparison_status"),
            "strict_acceptance": str((output / "strict_acceptance.json").resolve()),
        }, indent=2))
        return 0 if (
            manifest.get("status") == "pass"
            and manifest.get("comparison_status") == "pass"
        ) else 1
    if output.exists() and any(output.iterdir()):
        parser.error("--output-dir must not exist or must be empty; reuse is forbidden")
    if not args.runner.is_file():
        parser.error(f"runner does not exist: {args.runner}")
    output.mkdir(parents=True, exist_ok=True)
    manifest = workflow.materialize(
        args.generated_dir.resolve(), output, ["dd", "dg"]
    )
    manifest = patch_manifest(manifest, output, args.overlay.resolve())
    manifest = workflow.execute(manifest, args.runner.resolve(), output, None)
    expected_comparisons = {
        f"{branch}_{curve}_curve"
        for branch in ("dd", "dg") for curve in ("idvg", "idvd")
    }
    if set(manifest.get("comparisons", {})) != expected_comparisons:
        manifest["comparison_status"] = "incomplete"
        workflow.write_json(output / "workflow_manifest.json", manifest)
    else:
        apply_strict_acceptance(manifest, output)
    write_summary(manifest, output)
    print(json.dumps({
        "status": manifest.get("status"),
        "comparison_status": manifest.get("comparison_status"),
        "manifest": str((output / "workflow_manifest.json").resolve()),
        "summary": str((output / "continuous_baseline_summary.md").resolve()),
    }, indent=2))
    return 0 if (
        manifest.get("status") == "pass"
        and manifest.get("comparison_status") == "pass"
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
