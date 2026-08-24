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
import json
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
    (run_dir / "continuous_baseline_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--overlay", type=Path, default=DEFAULT_OVERLAY)
    args = parser.parse_args()
    output = args.output_dir.resolve()
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
