#!/usr/bin/env python3
"""Aggregate the prospective PN2D node-volume default-policy evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

try:
    from scripts.generate_pn2d_config import render_named_template
except ModuleNotFoundError:
    from generate_pn2d_config import render_named_template

try:
    from scripts.analyze_pn2d_avalanche_on_bv_parity import (
        KNEE_BIASES_V,
        adjacent_slopes,
        continuous_breakpoint,
        load_curve,
        slope_knee,
    )
except ModuleNotFoundError:
    from analyze_pn2d_avalanche_on_bv_parity import (
        KNEE_BIASES_V,
        adjacent_slopes,
        continuous_breakpoint,
        load_curve,
        slope_knee,
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def default_render_binding(
    controls: dict[str, Any], template_path: Path, contract: dict[str, Any]
) -> dict[str, Any]:
    manifest_record = controls.get("base_config_manifest") or {}
    candidate_record = controls.get("candidate_config") or {}
    manifest_path = Path(str(manifest_record.get("path", "")))
    config_path = Path(str(candidate_record.get("path", "")))
    gates: dict[str, bool] = {
        "template_hash_matches_contract": (
            template_path.is_file()
            and sha256(template_path)
            == contract["artifact_bindings"]["pn2d_bv_template_sha256"]
        ),
        "candidate_origin_is_default_render": (
            controls.get("candidate_origin") == "pn2d_bv_template_default_render"
        ),
        "driver_atomic_binding_passed": all(
            (controls.get("default_render_binding") or {}).get("gates", {}).values()
        ),
        "manifest_exists": manifest_path.is_file(),
        "config_exists": config_path.is_file(),
    }
    manifest: dict[str, Any] = {}
    if gates["manifest_exists"]:
        manifest = load_json(manifest_path)
        gates.update(
            {
                "manifest_hash_bound": (
                    sha256(manifest_path) == manifest_record.get("sha256")
                ),
                "manifest_has_no_profile_override": (
                    "avalanche_current_support_profile"
                    not in manifest.get("overrides", {})
                ),
                "manifest_version_qualified": (
                    int(manifest.get("template_version", 0))
                    >= int(contract["candidate"]["required_template_version_minimum"])
                ),
            }
        )
    if gates["manifest_exists"] and gates["config_exists"]:
        expected_config, expected_manifest = render_named_template(
            "pn2d_bv", manifest.get("overrides", {}), allow_absolute_paths=True
        )
        gates.update(
            {
                "manifest_replays_exactly": manifest == expected_manifest,
                "config_replays_exactly": load_json(config_path) == expected_config,
                "config_hash_bound": sha256(config_path) == candidate_record.get("sha256"),
            }
        )
    return {
        "passed": bool(gates) and all(gates.values()),
        "gates": gates,
        "manifest": (
            {"path": str(manifest_path.resolve()), "sha256": sha256(manifest_path)}
            if manifest_path.is_file()
            else None
        ),
        "config": (
            {"path": str(config_path.resolve()), "sha256": sha256(config_path)}
            if config_path.is_file()
            else None
        ),
    }


def scope_and_rollback_gates(
    bv_template_path: Path, iv_template_path: Path, contract: dict[str, Any]
) -> dict[str, Any]:
    default_bv, default_manifest = render_named_template("pn2d_bv")
    rollback_bv, rollback_manifest = render_named_template(
        "pn2d_bv",
        {"avalanche_current_support_profile": "legacy_cell_reconstructed"},
    )
    iv_config, _ = render_named_template("pn2d_iv")

    def observed_profile(config: dict[str, Any]) -> dict[str, Any]:
        impact = config["solver"]["impact_ionization"]
        return {
            "impact_ionization": {
                name: impact.get(name)
                for name in (
                    "current_approximation",
                    "source_mapping_mode",
                    "cell_reconstructed_midpoint_density",
                )
            },
            "mesh_geometry": config.get("mesh_geometry"),
        }

    gates = {
        "pn2d_bv_template_hash_matches_contract": (
            sha256(bv_template_path)
            == contract["artifact_bindings"]["pn2d_bv_template_sha256"]
        ),
        "pn2d_iv_template_hash_matches_contract": (
            sha256(iv_template_path)
            == contract["artifact_bindings"]["pn2d_iv_template_sha256"]
        ),
        "default_profile_is_atomic_candidate": (
            observed_profile(default_bv)
            == {
                "impact_ionization": contract["candidate"]["impact_ionization"],
                "mesh_geometry": contract["candidate"]["mesh_geometry"],
            }
        ),
        "default_render_has_no_profile_override": (
            "avalanche_current_support_profile"
            not in default_manifest.get("overrides", {})
        ),
        "rollback_profile_is_atomic": (
            observed_profile(rollback_bv)
            == {
                "impact_ionization": contract["rollback"]["impact_ionization"],
                "mesh_geometry": contract["rollback"]["mesh_geometry"],
            }
        ),
        "rollback_is_explicit_override": (
            rollback_manifest.get("overrides", {}).get(
                "avalanche_current_support_profile"
            )
            == "legacy_cell_reconstructed"
        ),
        "pn2d_iv_has_no_bv_mesh_policy": "mesh_geometry" not in iv_config,
        "pn2d_iv_impact_ionization_remains_off": (
            iv_config["solver"]["impact_ionization"].get("model") == "none"
        ),
    }
    return {"passed": all(gates.values()), "gates": gates}


def on_metrics(vela_csv: Path, sentaurus_csv: Path) -> dict[str, Any]:
    vela = {round(point.bias_V, 9): point for point in load_curve(vela_csv)}
    sentaurus = {
        round(point.bias_V, 9): point for point in load_curve(sentaurus_csv)
    }
    biases = sorted((set(vela) & set(sentaurus)) - {0.0}, reverse=True)
    knee = [bias for bias in KNEE_BIASES_V if bias in vela and bias in sentaurus]

    def errors(selected: list[float]) -> dict[str, float]:
        values = [
            abs(
                math.log10(
                    abs(
                        vela[bias].current_A_per_um
                        / sentaurus[bias].current_A_per_um
                    )
                )
            )
            for bias in selected
        ]
        return {
            "count": len(values),
            "rmse_dex": math.sqrt(sum(value * value for value in values) / len(values)),
            "maximum_dex": max(values),
        }

    ordered = sorted(set(vela) & set(sentaurus), reverse=True)
    nonmonotonic = [
        [ordered[index - 1], ordered[index]]
        for index in range(1, len(ordered))
        if abs(vela[ordered[index]].current_A_per_um)
        < abs(vela[ordered[index - 1]].current_A_per_um)
    ]
    vela_slopes = adjacent_slopes(vela, knee, "vela")
    sentaurus_slopes = adjacent_slopes(sentaurus, knee, "sentaurus")
    return {
        "all_nonzero": errors(biases),
        "knee": errors(knee),
        "V_break_V": {
            "vela": continuous_breakpoint(vela, knee, "vela"),
            "sentaurus": continuous_breakpoint(sentaurus, knee, "sentaurus"),
        },
        "V_slope_V": {
            "vela": slope_knee(vela_slopes),
            "sentaurus": slope_knee(sentaurus_slopes),
        },
        "nonmonotonic_intervals_V": nonmonotonic,
    }


def on_gate(metrics: dict[str, Any], gate: dict[str, Any]) -> dict[str, bool]:
    vela_slope = metrics["V_slope_V"]["vela"]
    sentaurus_slope = metrics["V_slope_V"]["sentaurus"]
    typed_slope = (vela_slope is None and sentaurus_slope is None) or (
        vela_slope is not None and sentaurus_slope is not None
    )
    return {
        "all_rmse": metrics["all_nonzero"]["rmse_dex"]
        <= gate["maximum_all_nonzero_log10_current_rmse_dex"],
        "all_maximum": metrics["all_nonzero"]["maximum_dex"]
        <= gate["maximum_all_nonzero_log10_current_error_dex"],
        "knee_rmse": metrics["knee"]["rmse_dex"]
        <= gate["maximum_knee_log10_current_rmse_dex"],
        "knee_maximum": metrics["knee"]["maximum_dex"]
        <= gate["maximum_knee_log10_current_error_dex"],
        "V_break": abs(
            metrics["V_break_V"]["vela"] - metrics["V_break_V"]["sentaurus"]
        )
        <= gate["maximum_V_break_absolute_error_V"],
        "V_slope_typed_outcome": typed_slope,
        "monotonicity": not metrics["nonmonotonic_intervals_V"],
    }


def release_gate(log_path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    text = log_path.read_text(encoding="utf-8-sig")
    summary = re.search(r"100% tests passed out of (\d+)", text)
    focused = {
        name: bool(re.search(re.escape(name) + r".*Passed", text))
        for name in contract["release_gate"]["required_focused_tests"]
    }
    return {
        "test_count": int(summary.group(1)) if summary else 0,
        "zero_failed": summary is not None,
        "focused_tests": focused,
        "passed": summary is not None and all(focused.values()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--iv-template", type=Path, required=True)
    parser.add_argument("--m0-root", type=Path, required=True)
    parser.add_argument("--m0-sentaurus-root", type=Path, required=True)
    parser.add_argument("--m2-root", type=Path, required=True)
    parser.add_argument("--m2-sentaurus-root", type=Path, required=True)
    parser.add_argument("--forward-report", type=Path, required=True)
    parser.add_argument("--release-log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract_path = args.contract.resolve()
    contract = load_json(contract_path)
    scope_guards = scope_and_rollback_gates(
        args.template.resolve(), args.iv_template.resolve(), contract
    )
    grids: dict[str, Any] = {}
    for name, root_arg, sentaurus_arg in (
        ("M0", args.m0_root, args.m0_sentaurus_root),
        ("M2", args.m2_root, args.m2_sentaurus_root),
    ):
        root = root_arg.resolve()
        sentaurus = sentaurus_arg.resolve()
        controls_path = root / "gate_report.json"
        controls = load_json(controls_path)
        binding = default_render_binding(
            controls, args.template.resolve(), contract
        )
        metrics = on_metrics(
            root / "avalanche_on" / "run-a" / "avalanche_on" / "iv.csv",
            sentaurus / "avalanche_on" / "normalized" / "aggregate.csv",
        )
        gates = on_gate(metrics, contract["avalanche_on_gate"])
        grids[name] = {
            "controls_path": str(controls_path),
            "controls_sha256": sha256(controls_path),
            "controls_passed": controls.get("status") == "passed",
            "default_render_binding": binding,
            "avalanche_on_metrics": metrics,
            "avalanche_on_gates": gates,
            "passed": (
                controls.get("status") == "passed"
                and binding["passed"]
                and all(gates.values())
            ),
        }
    forward_path = args.forward_report.resolve()
    forward = load_json(forward_path)
    release_path = args.release_log.resolve()
    release = release_gate(release_path, contract)
    passed = (
        all(grid["passed"] for grid in grids.values())
        and scope_guards["passed"]
        and forward.get("status") == "passed"
        and release["passed"]
    )
    report = {
        "schema": "vela.pn2d_node_volume_policy_default_acceptance_review.v1",
        "status": "passed" if passed else "failed",
        "outcome": (
            "ready_for_independent_default_policy_reviews"
            if passed
            else "default_policy_acceptance_gate_failed"
        ),
        "contract": {"path": str(contract_path), "sha256": sha256(contract_path)},
        "scope_and_rollback": scope_guards,
        "grids": grids,
        "forward_iv": {
            "path": str(forward_path),
            "sha256": sha256(forward_path),
            "passed": forward.get("status") == "passed",
            "anchors": forward.get("anchors"),
        },
        "release": {
            "path": str(release_path),
            "sha256": sha256(release_path),
            **release,
        },
        "authorization": {
            "production_default_change_authorized": False,
            "independent_scientific_review_required": True,
            "independent_code_review_required": True,
        },
    }
    args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
    args.output.resolve().write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
