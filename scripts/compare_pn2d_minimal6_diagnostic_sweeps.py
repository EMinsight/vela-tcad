#!/usr/bin/env python3
"""Compare exact PN2D Minimal6 diagnostic-sweep checkpoints without interpolation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, UnidentifiedImageError

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.pn2d_minimal6_diagnostics.schemas import (
    DISCLAIMER,
    validate_bv_comparison_v1,
    validate_formula_difference_v1,
    validate_sweep_manifest_v1,
)
from scripts.run_pn2d_minimal6_diagnostic_sweep import (
    validate_sweep_manifest,
    BRANCH_THRESHOLD_VERSION,
    GEOMETRIC_ZERO_SOURCE_INTEGRAL,
    classify_branch,
)

SCHEMA = "vela.pn2d_minimal6_bv_comparison.v1"
EPSILON = 1.0e-12
CONTACT_TOLERANCE_RELATIVE = 1.0e-9
CONTACT_TOLERANCE_FLOOR_A_PER_UM = 1.0e-18
SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")
OBSERVABLES = (
    "anode_current_A_per_um", "cathode_current_A_per_um", "max_field_V_per_m",
    "native_source_integral_s_inv_per_cm", "reconstructed_source_integral_s_inv_per_cm",
)
FIGURE_SCHEMA = "vela.pn2d_minimal6_figure_contract.v1"
FIGURE_NAMES = (
    "terminal_current.png",
    "one_volt_growth.png",
    "maximum_field.png",
    "source_integrals.png",
    "topology.png",
)
GENERATED_ARTIFACT_NAMES = {
    "sweep_comparison.csv", "sweep_comparison.md", *FIGURE_NAMES,
}
INPUT_ARTIFACT_NAMES = {
    "vela_manifest", "sentaurus_manifest", "fixed_state_report",
}
EXPECTED_ARTIFACT_HASH_NAMES = GENERATED_ARTIFACT_NAMES | {
    f"input:{name}" for name in INPUT_ARTIFACT_NAMES
}
POST_GENERATION_FIELDS = {"figure_contract", "artifact_hashes", "input_artifacts"}
FIGURE_WIDTH_PX = 900
FIGURE_HEIGHT_PX = 504
FIGURE_DPI = 120
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
FIGURE_METADATA = {
    "DiagnosticDisclaimer": DISCLAIMER,
    "SolverTermination": "Every recorded solver failure transition is explicitly marked.",
    "BVExtrapolation": "No physical breakdown voltage (BV) is extrapolated.",
}


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_sha(payload: Any) -> str:
    return _sha_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _finite(value: Any) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(float(value))

def _valid_sha256(value: Any) -> bool:
    return isinstance(value, str) and SHA256.fullmatch(value) is not None


def ratio_record(numerator: float | None, denominator: float | None) -> dict[str, Any]:
    """Return a typed ratio; zero and unavailable values are never coerced."""
    if isinstance(numerator, bool) or isinstance(denominator, bool):
        raise ValueError("boolean ratio inputs are forbidden")
    if numerator is None or denominator is None or not _finite(numerator) or not _finite(denominator):
        return {"classification": "unavailable", "value": None}
    if denominator == 0.0:
        return {"classification": "zero_denominator", "value": None}
    if numerator == 0.0:
        return {"classification": "zero_numerator", "value": 0.0}
    return {"classification": "available", "value": float(numerator / denominator)}


def _accepted(manifest: dict[str, Any], solver: str) -> list[dict[str, Any]]:
    rows = manifest.get("accepted_checkpoints", [])
    if not isinstance(rows, list):
        raise ValueError(f"{solver} manifest has invalid accepted_checkpoints")
    checked: list[dict[str, Any]] = []
    seen: set[tuple[str, float]] = set()
    for row in rows:
        if "quantity_ledger_result" in row:
            raise ValueError(
                f"{solver} checkpoint contains unverified embedded scientific summary: quantity_ledger_result")
        if row.get("solver") != solver or row.get("status") != "accepted":
            raise ValueError(f"{solver} manifest has invalid accepted checkpoint")
        actual, target = row.get("actual_bias_V"), row.get("target_bias_V")
        if not _finite(actual) or not _finite(target) or abs(float(actual) - float(target)) > EPSILON:
            raise ValueError(f"{solver} checkpoint is not an exact target bias")
        state_path, state_sha256 = row.get("state_path"), row.get("state_sha256")
        if not isinstance(state_path, str) or not state_path or not _valid_sha256(state_sha256):
            raise ValueError(f"{solver} checkpoint lacks hash-addressed state identity")
        topology = row.get("topology")
        if not isinstance(topology, str) or not topology:
            raise ValueError(f"{solver} checkpoint lacks topology")
        key = (topology, float(target))
        if key in seen:
            raise ValueError(f"{solver} has duplicate exact checkpoint {key}")
        seen.add(key)
        observables = row.get("observables")
        if not isinstance(observables, dict) or any(not _finite(observables.get(name)) for name in OBSERVABLES):
            raise ValueError(f"{solver} checkpoint lacks finite observables")
        if row.get("branch_classification") != "unidentified":
            raise ValueError(f"{solver} checkpoint has noncanonical side-only branch evidence")
        threshold = manifest.get("branch_threshold_version")
        if threshold != BRANCH_THRESHOLD_VERSION or row.get("branch_threshold_version") != threshold:
            raise ValueError(f"{solver} checkpoint lacks canonical branch-threshold provenance")
        convergence = row.get("convergence_metadata")
        if not isinstance(convergence, dict) or not convergence:
            raise ValueError(f"{solver} checkpoint lacks convergence metadata")
        checked.append(row)
    return checked


def _index(rows: list[dict[str, Any]]) -> dict[tuple[str, float], dict[str, Any]]:
    return {(str(row["topology"]), float(row["target_bias_V"])): row for row in rows}


def _current(row: dict[str, Any]) -> float:
    return float(row["observables"]["anode_current_A_per_um"])


def _ratio_rows(numerator: dict[str, Any], denominator: dict[str, Any], observable: str, *, absolute: bool = False) -> dict[str, Any]:
    left = float(numerator["observables"][observable])
    right = float(denominator["observables"][observable])
    return ratio_record(abs(left) if absolute else left, abs(right) if absolute else right)


def _source_geometric_zero_pair(
    numerator: dict[str, Any], denominator: dict[str, Any], observable: str,
) -> bool:
    return all(
        abs(float(row["observables"][observable])) <= GEOMETRIC_ZERO_SOURCE_INTEGRAL
        for row in (numerator, denominator)
    )


def _branch_geometric_zero_pair(*rows: dict[str, Any]) -> bool:
    source_fields = (
        "native_source_integral_s_inv_per_cm",
        "reconstructed_source_integral_s_inv_per_cm",
    )
    return all(
        abs(float(row["observables"][field])) <= GEOMETRIC_ZERO_SOURCE_INTEGRAL
        for row in rows for field in source_fields
    )


def _comparison_ratio(
    numerator: dict[str, Any],
    denominator: dict[str, Any],
    observable: str,
    *,
    absolute: bool = False,
    geometric_zero: bool = False,
) -> dict[str, Any]:
    if geometric_zero:
        return {"classification": "geometric_zero", "value": None}
    return _ratio_rows(numerator, denominator, observable, absolute=absolute)


def _contact_unavailable() -> dict[str, Any]:
    return {
        "classification": "unavailable", "value": None,
        "reason": "contact_current_not_conserved",
    }


def _is_contact_conserved(row: dict[str, Any]) -> bool:
    return _contact_conservation(row)["classification"] == "conserved"


def _one_volt_growth(index: dict[tuple[str, float], dict[str, Any]], topology: str, bias: float) -> dict[str, Any]:
    current, next_row = index.get((topology, bias)), index.get((topology, bias - 1.0))
    if current is None or next_row is None:
        return {"classification": "unavailable", "value": None, "reason": "next exact one-volt checkpoint unavailable"}
    if not _is_contact_conserved(current) or not _is_contact_conserved(next_row):
        return _contact_unavailable()
    if _current(current) == 0.0 or _current(next_row) == 0.0:
        return {"classification": "zero_current", "value": None}
    return ratio_record(abs(_current(next_row)), abs(_current(current)))


def _failures(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    rows = manifest.get("failed_transitions", [])
    if not isinstance(rows, list):
        raise ValueError("failed_transitions must be a list")
    for row in rows:
        if row.get("status") != "rejected" or row.get("observables") is not None:
            raise ValueError("failure transition has fabricated observables")
        if row.get("branch_classification") != "unidentified":
            raise ValueError("failure transition lacks unidentified branch classification")
        if row.get("branch_threshold_version") != BRANCH_THRESHOLD_VERSION:
            raise ValueError("failure transition lacks canonical branch-threshold provenance")
        convergence = row.get("convergence_metadata")
        if not isinstance(convergence, dict) or not convergence:
            raise ValueError("failure transition lacks convergence metadata")
        if not isinstance(row.get("incomplete_reason"), str) or not row["incomplete_reason"]:
            raise ValueError("failure transition lacks incomplete reason")
    return rows


def _gap_closure(ratios: dict[str, dict[str, Any]]) -> dict[str, Any]:
    gaps: list[dict[str, Any]] = []
    eligible = 0
    for quantity, ratio in ratios.items():
        value = ratio.get("value")
        is_eligible = (
            ratio.get("classification") == "available"
            and _finite(value) and float(value) > 0.0
        )
        if is_eligible:
            eligible += 1
        gaps.append({
            "quantity": quantity,
            "classification": str(ratio.get("classification", "unavailable")),
            "decomposition_status": "unidentifiable" if is_eligible else "not_applicable",
            "log_gap_dex": math.log10(abs(float(value))) if is_eligible else None,
            "named_contributions": [],
            "residual": {
                "name": "cross_solver_semantics_residual",
                "classification": "unidentifiable",
                "value_dex": None,
            },
            "closure_error_dex": None,
        })
    return {
        "status": "unidentifiable" if eligible else "not_applicable",
        "tolerance_dex": 1.0e-10,
        "gaps": gaps,
    }


def _validated_deck_rows(manifest: dict[str, Any], solver: str) -> list[dict[str, Any]]:
    deck_key = "segments" if solver == "vela" else "sentaurus_segments"
    rows = manifest.get(deck_key)
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{solver} manifest requires a nonempty {deck_key} deck set")
    checked: list[dict[str, Any]] = []
    identities: set[str] = set()
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"{solver} {deck_key}[{index}] must be a mapping")
        deck = row.get("deck")
        digest_value = row.get("deck_sha256")
        if not isinstance(deck, str) or not deck or not _valid_sha256(digest_value):
            raise ValueError(f"{solver} {deck_key}[{index}] lacks a hash-addressed deck")
        identity = str(Path(deck))
        if identity in identities:
            raise ValueError(f"{solver} {deck_key} contains duplicate deck artifact identity")
        identities.add(identity)
        checked.append(row)
    return checked


def _validate_manifest_configuration(manifest: dict[str, Any], solver: str) -> dict[str, Any]:
    template = manifest.get("template")
    if not isinstance(template, dict) or not _valid_sha256(template.get("sha256")):
        raise ValueError(f"{solver} manifest lacks a hash-addressed template")
    topology_hashes = manifest.get("topology_input_sha256")
    if not isinstance(topology_hashes, dict) or set(topology_hashes) != {"sketch", "mirror"}:
        raise ValueError(f"{solver} manifest requires sketch/mirror topology hashes")
    for topology, entries in topology_hashes.items():
        if not isinstance(entries, dict) or not entries:
            raise ValueError(f"{solver} {topology} topology hashes must be non-empty")
        if any(not isinstance(name, str) or not name or not _valid_sha256(value) for name, value in entries.items()):
            raise ValueError(f"{solver} {topology} topology hashes are invalid")
    deck_rows = _validated_deck_rows(manifest, solver)
    config = {
        "template": template,
        "topology_input_sha256": topology_hashes,
        "deck_sha256": sorted(row["deck_sha256"] for row in deck_rows),
        "sweep_manifest_sha256": _canonical_sha(manifest),
    }
    return {**config, "configuration_sha256": _canonical_sha(config)}


def _contact_conservation(row: dict[str, Any]) -> dict[str, Any]:
    observables = row["observables"]
    anode = float(observables["anode_current_A_per_um"])
    cathode = float(observables["cathode_current_A_per_um"])
    residual = anode + cathode
    scale = max(abs(anode), abs(cathode))
    tolerance = max(CONTACT_TOLERANCE_FLOOR_A_PER_UM, scale * CONTACT_TOLERANCE_RELATIVE)
    return {
        "anode_current_A_per_um": anode,
        "cathode_current_A_per_um": cathode,
        "signed_residual_A_per_um": residual,
        "tolerance_A_per_um": tolerance,
        "tolerance_formula": "max(1e-18 A/um, 1e-9 * max(|Anode|, |Cathode|))",
        "classification": "conserved" if abs(residual) <= tolerance else "not_conserved",
    }


def _checked_bound_artifact(
    package_root: Path,
    path_value: Any,
    expected_sha256: Any,
    label: str,
    *,
    allow_absolute: bool,
) -> Path:
    if not isinstance(path_value, str) or not path_value:
        raise ValueError(f"{label} lacks an artifact path")
    declared = Path(path_value)
    if declared.is_absolute():
        if not allow_absolute:
            raise ValueError(f"{label} path must be relative to the manifest root")
        resolved = declared.resolve()
    else:
        resolved = (package_root / declared).resolve()
        if resolved != package_root and package_root not in resolved.parents:
            raise ValueError(f"{label} relative path escapes the manifest root")
    if not resolved.is_file():
        raise ValueError(f"{label} artifact is missing: {resolved}")
    if not _valid_sha256(expected_sha256) or _sha_bytes(resolved.read_bytes()) != expected_sha256:
        raise ValueError(f"{label} artifact is hash-tampered: {resolved}")
    return resolved


def _validate_bound_sweep_input(
    manifest: dict[str, Any], manifest_path: Path, solver: str,
) -> None:
    validate_sweep_manifest_v1(manifest)
    path = Path(manifest_path).resolve()
    if not path.is_file():
        raise ValueError(f"{solver} sweep manifest input artifact is missing")
    root = path.parent
    template = manifest.get("template")
    if not isinstance(template, dict):
        raise ValueError(f"{solver} manifest lacks template provenance")
    _checked_bound_artifact(
        root, template.get("path"), template.get("sha256"),
        f"{solver} template", allow_absolute=True,
    )
    topology_hashes = manifest.get("topology_input_sha256")
    if not isinstance(topology_hashes, dict) or set(topology_hashes) != {"sketch", "mirror"}:
        raise ValueError(f"{solver} manifest requires sketch/mirror topology inputs")
    for topology in ("sketch", "mirror"):
        entries = topology_hashes[topology]
        if not isinstance(entries, dict) or not entries:
            raise ValueError(f"{solver} {topology} topology inputs must be nonempty")
        for name, expected_sha256 in entries.items():
            if (not isinstance(name, str) or not name or Path(name).name != name):
                raise ValueError(f"{solver} {topology} has an invalid declared input name")
            _checked_bound_artifact(
                root, str(Path("inputs") / topology / name), expected_sha256,
                f"{solver} {topology} input {name}", allow_absolute=False,
            )
    for index, row in enumerate(_validated_deck_rows(manifest, solver)):
        _checked_bound_artifact(
            root, row["deck"], row["deck_sha256"],
            f"{solver} deck row {index}", allow_absolute=False,
        )
    accepted = manifest.get("accepted_checkpoints")
    if not isinstance(accepted, list):
        raise ValueError(f"{solver} accepted checkpoints must be an array")
    for index, row in enumerate(accepted):
        if not isinstance(row, dict):
            raise ValueError(f"{solver} accepted checkpoint {index} must be a mapping")
        _checked_bound_artifact(
            root, row.get("state_path"), row.get("state_sha256"),
            f"{solver} accepted state {index}", allow_absolute=False,
        )
    if manifest.get("targets_V") == [float(-index) for index in range(21)]:
        validate_sweep_manifest(manifest, package_root=root)

def _fixed_state_recheck(common: dict[tuple[str, float], tuple[dict[str, Any], dict[str, Any]]], fixed: dict[str, Any]) -> list[dict[str, Any]]:
    dominance_rules = fixed.get("dominance_rules")
    if isinstance(dominance_rules, dict):
        original = dominance_rules.get("status", "unavailable")
        fixed_dominant = dominance_rules.get("dominant_factor")
    else:
        original = fixed.get("root_cause_status", "unavailable")
        fixed_dominant = fixed.get("dominant_factor")
    rows: list[dict[str, Any]] = []
    for bias in (0.0, -12.0, -19.0):
        topologies = [topology for topology in ("sketch", "mirror") if (topology, bias) in common]
        states: list[dict[str, Any]] = []
        for topology in topologies:
            vela, sentaurus = common[(topology, bias)]
            for solver, checkpoint in (("vela", vela), ("sentaurus", sentaurus)):
                states.append({
                    "solver": solver,
                    "topology": topology,
                    "bias_V": bias,
                    "state_path": checkpoint.get("state_path"),
                    "state_sha256": checkpoint.get("state_sha256"),
                    "state_binding_status": "manifest_hash_addressed",
                })
        row: dict[str, Any] = {
            "bias_V": bias,
            "status": "unidentifiable",
            "ranking_status": "unidentifiable",
            "fixed_state_status": original,
            "fixed_state_dominant_factor": fixed_dominant,
            "topologies": topologies,
            "self_consistent_states": states,
        }
        if states:
            row.update({
                "reason_code": "missing_verified_nonlinear_ledger_input_bundle",
                "reason": "verified nonlinear ledger-input bundle is unavailable for the hash-addressed self-consistent checkpoints",
                "missing_inputs": ["verified_nonlinear_ledger_input_bundle"],
                "recheck_basis": "hash_addressed_self_consistent_states_without_verified_nonlinear_ledger_bundle",
            })
        else:
            row.update({
                "reason_code": "no_common_self_consistent_state",
                "reason": "no common exact self-consistent checkpoints are available at this fixed-state bias",
                "missing_inputs": [
                    "common_exact_self_consistent_states",
                    "verified_nonlinear_ledger_input_bundle",
                ],
                "recheck_basis": "no_common_exact_self_consistent_states",
            })
        rows.append(row)
    return rows


def compare_sweeps(vela_manifest: dict[str, Any], sentaurus_manifest: dict[str, Any], *, fixed_state_report: dict[str, Any]) -> dict[str, Any]:
    """Produce a comparison object from accepted, exact-bias rows only."""
    validate_sweep_manifest_v1(vela_manifest)
    validate_sweep_manifest_v1(sentaurus_manifest)
    vela_threshold = vela_manifest.get("branch_threshold_version")
    sentaurus_threshold = sentaurus_manifest.get("branch_threshold_version")
    if not isinstance(vela_threshold, str) or not vela_threshold or vela_threshold != sentaurus_threshold:
        raise ValueError("sweep manifests must declare the same non-empty branch threshold version")
    branch_threshold_version = vela_threshold
    vela_configuration = _validate_manifest_configuration(vela_manifest, "vela")
    sentaurus_configuration = _validate_manifest_configuration(sentaurus_manifest, "sentaurus")
    vela_rows, sentaurus_rows = _accepted(vela_manifest, "vela"), _accepted(sentaurus_manifest, "sentaurus")
    vela_index, sentaurus_index = _index(vela_rows), _index(sentaurus_rows)
    keys = sorted(set(vela_index) & set(sentaurus_index), key=lambda key: (key[0], -key[1]))
    common = {key: (vela_index[key], sentaurus_index[key]) for key in keys}
    checkpoints: list[dict[str, Any]] = []
    for topology, bias in keys:
        vela, sentaurus = common[(topology, bias)]
        branch_geometric_zero = _branch_geometric_zero_pair(vela, sentaurus)
        sentaurus_current = _current(sentaurus)
        vela_current = _current(vela)
        zero_current = sentaurus_current == 0.0 or vela_current == 0.0
        absolute_ratio = None if sentaurus_current == 0.0 else abs(vela_current / sentaurus_current)
        branch_ratio_evidence = {
            "vela_anode_current_A_per_um": vela_current,
            "sentaurus_anode_current_A_per_um": sentaurus_current,
            "absolute_vela_over_sentaurus": absolute_ratio,
            "geometric_zero": branch_geometric_zero,
            "threshold_version": branch_threshold_version,
        }
        contact_evidence = {
            "vela": _contact_conservation(vela),
            "sentaurus": _contact_conservation(sentaurus),
        }
        contacts_conserved = all(
            evidence["classification"] == "conserved"
            for evidence in contact_evidence.values()
        )
        ratios = {
            "terminal_current": (
                _comparison_ratio(vela, sentaurus, "anode_current_A_per_um", absolute=True)
                if contacts_conserved else _contact_unavailable()
            ),
            "maximum_field": _comparison_ratio(vela, sentaurus, "max_field_V_per_m"),
            "native_source": _comparison_ratio(
                vela, sentaurus, "native_source_integral_s_inv_per_cm",
                geometric_zero=_source_geometric_zero_pair(
                    vela, sentaurus, "native_source_integral_s_inv_per_cm"
                ),
            ),
            "reconstructed_source": _comparison_ratio(
                vela, sentaurus, "reconstructed_source_integral_s_inv_per_cm",
                geometric_zero=_source_geometric_zero_pair(
                    vela, sentaurus, "reconstructed_source_integral_s_inv_per_cm"
                ),
            ),
        }
        if not contacts_conserved:
            sign_alignment = _contact_unavailable()
        elif zero_current:
            sign_alignment = {"classification": "zero_current", "value": None}
        else:
            sign_alignment = {
                "classification": "available",
                "value": "aligned" if math.copysign(1.0, vela_current) == math.copysign(1.0, sentaurus_current) else "opposed",
            }
        branch_classification = (
            "unidentified"
            if not contacts_conserved
            else classify_branch(
                sentaurus_current, vela_current,
                geometric_zero=branch_geometric_zero or zero_current,
            )
        )
        checkpoints.append({
            "topology": topology, "bias_V": bias, "classification": "common_exact", "vela": vela, "sentaurus": sentaurus,
            "branch_classification": branch_classification,
            "branch_threshold_version": branch_threshold_version,
            "branch_ratio_evidence": branch_ratio_evidence,
            "contact_current_conservation": {
                "unit": "A/um", **contact_evidence,
            },
            "terminal_current_sign_alignment": sign_alignment,
            "terminal_current_ratio": ratios["terminal_current"],
            "maximum_field_ratio": ratios["maximum_field"],
            "native_source_ratio": ratios["native_source"],
            "reconstructed_source_ratio": ratios["reconstructed_source"],
            "vela_one_volt_current_growth": _one_volt_growth(vela_index, topology, bias),
            "sentaurus_one_volt_current_growth": _one_volt_growth(sentaurus_index, topology, bias),
            "gap_closure": _gap_closure(ratios),
        })
    side_only: list[dict[str, Any]] = []
    missing_tails: list[dict[str, Any]] = []
    for solver, index, other in (("vela", vela_index, sentaurus_index), ("sentaurus", sentaurus_index, vela_index)):
        by_topology: dict[str, list[float]] = {}
        for (topology, bias), row in index.items():
            if (topology, bias) not in other:
                side_only.append({"solver": solver, "topology": topology, "bias_V": bias, "classification": "side_only", "checkpoint": row})
                by_topology.setdefault(topology, []).append(bias)
        for topology, biases in by_topology.items():
            missing_tails.append({"solver": solver, "topology": topology, "biases_V": sorted(biases, reverse=True), "reason": "no accepted exact checkpoint from the other solver"})
    topology_sensitivity: list[dict[str, Any]] = []
    for solver, index in (("vela", vela_index), ("sentaurus", sentaurus_index)):
        biases = sorted({bias for topology, bias in index if topology == "sketch" and ("mirror", bias) in index}, reverse=True)
        for bias in biases:
            sketch, mirror = index[("sketch", bias)], index[("mirror", bias)]
            topology_sensitivity.append({"solver": solver, "bias_V": bias,
                "terminal_current_sketch_over_mirror": (
                    ratio_record(abs(_current(sketch)), abs(_current(mirror)))
                    if _is_contact_conserved(sketch) and _is_contact_conserved(mirror)
                    else _contact_unavailable()
                ),
                "maximum_field_sketch_over_mirror": _ratio_rows(sketch, mirror, "max_field_V_per_m"),
                "native_source_sketch_over_mirror": _comparison_ratio(
                    sketch, mirror, "native_source_integral_s_inv_per_cm",
                    geometric_zero=_source_geometric_zero_pair(
                        sketch, mirror, "native_source_integral_s_inv_per_cm"),
                )})
    failures = _failures(vela_manifest) + _failures(sentaurus_manifest)
    deepest = min((bias for _, bias in keys), default=None)
    eligible_gap_count = sum(
        1 for checkpoint in checkpoints for gap in checkpoint["gap_closure"]["gaps"]
        if gap["log_gap_dex"] is not None
    )
    report = {
        "schema": SCHEMA, "diagnostic_disclaimer": DISCLAIMER, "interpolation": "forbidden",
        "comparison_status": "available" if keys else "stopped_with_evidence",
        "validation_failure": None if keys else {"code": "no_exact_common_checkpoint", "message": "no accepted identity-verified exact checkpoint is common to both solvers"},
        "branch_threshold_version": branch_threshold_version,
        "solver_configurations": {
            "vela": vela_configuration,
            "sentaurus": sentaurus_configuration,
        },
        "accepted_transitions": {"vela": vela_rows, "sentaurus": sentaurus_rows},
        "failed_transitions": failures, "failure_transitions": failures,
        "checkpoints": checkpoints, "records": checkpoints, "terminal_currents": checkpoints, "maximum_fields": checkpoints, "source_integrals": checkpoints,
        "convergence_metadata": {"vela_accepted": len(vela_rows), "sentaurus_accepted": len(sentaurus_rows), "common_exact": len(checkpoints)},
        "curve_artifact_hashes": {"vela_manifest": _canonical_sha(vela_manifest), "sentaurus_manifest": _canonical_sha(sentaurus_manifest)},
        "deepest_common_bias_V": {"classification": "available", "value": deepest} if deepest is not None else {"classification": "unavailable", "value": None, "reason": "no accepted exact common checkpoint"},
        "missing_tails": sorted(missing_tails, key=lambda row: (row["solver"], row["topology"])),
        "side_only_checkpoints": sorted(side_only, key=lambda row: (row["solver"], row["topology"], -row["bias_V"])),
        "topology_sensitivity": topology_sensitivity,
        "fixed_state_recheck": _fixed_state_recheck(common, fixed_state_report),
        "artifact_hashes": {},
        "input_artifacts": {},
        "closure": {
            "status": "unidentifiable" if eligible_gap_count else "not_applicable",
            "eligible_gaps": eligible_gap_count,
            "decomposed_gaps": 0,
            "unidentifiable_gaps": eligible_gap_count,
            "rule": "observed positive log gaps are retained without fabricated decomposition",
        },
    }
    return report


def _finish(fig: plt.Figure, ax: plt.Axes, title: str, ylabel: str) -> None:
    ax.set_title(title); ax.set_ylabel(ylabel); ax.set_xlabel("exact applied bias (V)"); ax.grid(True, alpha=0.25)
    fig.text(0.01, 0.01, DISCLAIMER + "; solver termination is marked; no BV extrapolation", fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1))


def _side_plot(ax: plt.Axes, report: dict[str, Any], observable: str, title: str, ylabel: str, absolute: bool = False) -> None:
    any_data = False
    for solver in ("vela", "sentaurus"):
        rows = report["accepted_transitions"][solver]
        for topology in ("sketch", "mirror"):
            selected = sorted((row for row in rows if row["topology"] == topology), key=lambda row: row["target_bias_V"], reverse=True)
            if selected:
                x = [row["target_bias_V"] for row in selected]; y = [float(row["observables"][observable]) for row in selected]
                ax.plot(x, [abs(value) for value in y] if absolute else y, marker="o", label=f"{solver} {topology}"); any_data = True
    _mark_failure_transitions(ax, report)
    if any_data: ax.legend(loc="upper left")
    else: ax.text(0.5, 0.5, "No accepted checkpoint", ha="center", va="center", transform=ax.transAxes)
    _finish(ax.figure, ax, title, ylabel)


def _mark_failure_transitions(ax: plt.Axes, report: dict[str, Any]) -> None:
    """Mark each recorded solver termination on one diagnostic figure."""
    for failure in report["failure_transitions"]:
        ax.axvline(float(failure["target_bias_V"]), color="tab:red", linestyle="--", alpha=0.55)



def _failure_marker_identities(report: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "solver": row["solver"],
            "topology": row["topology"],
            "start_bias_V": row["start_bias_V"],
            "target_bias_V": row["target_bias_V"],
        }
        for row in report["failure_transitions"]
    ]


def _accepted_series_identities(report: dict[str, Any], quantity: str) -> list[dict[str, Any]]:
    identities: list[dict[str, Any]] = []
    for solver in ("vela", "sentaurus"):
        rows = report["accepted_transitions"][solver]
        for topology in ("sketch", "mirror"):
            if any(row["topology"] == topology for row in rows):
                identities.append({"solver": solver, "topology": topology, "quantity": quantity})
    return identities


def _save_figure(
    fig: plt.Figure,
    path: Path,
    *,
    series_identities: list[dict[str, Any]],
    failure_markers: list[dict[str, Any]],
) -> dict[str, Any]:
    metadata = {
        **FIGURE_METADATA,
        "SeriesIdentities": json.dumps(series_identities, sort_keys=True, separators=(",", ":")),
        "FailureTransitionMarkers": json.dumps(failure_markers, sort_keys=True, separators=(",", ":")),
    }
    fig.savefig(path, dpi=FIGURE_DPI, metadata=metadata)
    plt.close(fig)
    payload = path.read_bytes()
    return {
        "sha256": _sha_bytes(payload),
        "width_px": FIGURE_WIDTH_PX,
        "height_px": FIGURE_HEIGHT_PX,
        "metadata": dict(FIGURE_METADATA),
        "series_identities": series_identities,
        "failure_transition_markers": failure_markers,
    }


def _render_figures(out_dir: Path, report: dict[str, Any]) -> dict[str, Any]:
    entries: dict[str, dict[str, Any]] = {}
    failure_markers = _failure_marker_identities(report)
    for stem, observable, title, unit, absolute, quantity in (
        ("terminal_current", "anode_current_A_per_um", "Terminal current at accepted exact checkpoints", "A/um", False, "terminal_current"),
        ("maximum_field", "max_field_V_per_m", "Maximum electric field at accepted exact checkpoints", "V/m", False, "maximum_field"),

    ):
        fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=FIGURE_DPI)
        _side_plot(ax, report, observable, title, unit, absolute)
        name = f"{stem}.png"
        entries[name] = _save_figure(
            fig,
            out_dir / name,
            series_identities=_accepted_series_identities(report, quantity),
            failure_markers=failure_markers,
        )

    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=FIGURE_DPI)
    source_series: list[dict[str, Any]] = []
    source_data = False
    for solver in ("vela", "sentaurus"):
        rows = report["accepted_transitions"][solver]
        for topology in ("sketch", "mirror"):
            selected = sorted((row for row in rows if row["topology"] == topology), key=lambda row: row["target_bias_V"], reverse=True)
            for observable, quantity, label in (
                ("native_source_integral_s_inv_per_cm", "native_source", "native source"),
                ("reconstructed_source_integral_s_inv_per_cm", "reconstructed_source", "reconstructed source"),
            ):
                if selected:
                    ax.plot(
                        [row["target_bias_V"] for row in selected],
                        [abs(float(row["observables"][observable])) for row in selected],
                        marker="o",
                        label=f"{solver} {topology} {label}",
                    )
                    source_series.append({"solver": solver, "topology": topology, "quantity": quantity})
                    source_data = True
    _mark_failure_transitions(ax, report)
    if source_data:
        ax.legend(loc="upper left")
    else:
        ax.text(0.5, 0.5, "No accepted checkpoint", ha="center", va="center", transform=ax.transAxes)
    _finish(fig, ax, "Native and reconstructed avalanche sources", "s^-1 per 1 cm depth")
    entries["source_integrals.png"] = _save_figure(
        fig,
        out_dir / "source_integrals.png",
        series_identities=source_series,
        failure_markers=failure_markers,
    )

    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=FIGURE_DPI)
    growth_series: list[dict[str, Any]] = []
    for solver, field, marker, label in (
        ("vela", "vela_one_volt_current_growth", "o", "Vela"),
        ("sentaurus", "sentaurus_one_volt_current_growth", "s", "Sentaurus"),
    ):
        selected = [row for row in report["checkpoints"] if row[field]["classification"] == "available"]
        if selected:
            ax.plot(
                [row["bias_V"] for row in selected],
                [row[field]["value"] for row in selected],
                marker=marker,
                label=label,
            )
            growth_series.append({"solver": solver, "quantity": "one_volt_growth"})
    if growth_series:
        ax.legend(loc="upper left")
    else:
        ax.text(0.5, 0.5, "No exact one-volt pair", ha="center", va="center", transform=ax.transAxes)
    _mark_failure_transitions(ax, report)
    _finish(fig, ax, "One-volt terminal-current growth", "growth ratio")
    entries["one_volt_growth.png"] = _save_figure(
        fig,
        out_dir / "one_volt_growth.png",
        series_identities=growth_series,
        failure_markers=failure_markers,
    )

    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=FIGURE_DPI)
    topology = [row for row in report["topology_sensitivity"] if row["terminal_current_sketch_over_mirror"]["classification"] == "available"]
    topology_series: list[dict[str, Any]] = []
    if topology:
        for solver in ("vela", "sentaurus"):
            selected = [row for row in topology if row["solver"] == solver]
            if selected:
                ax.plot([row["bias_V"] for row in selected], [row["terminal_current_sketch_over_mirror"]["value"] for row in selected], marker="o", label=solver)
                topology_series.append({"solver": solver, "quantity": "terminal_current_sketch_over_mirror"})
        ax.legend(loc="upper left")
    else:
        ax.text(0.5, 0.5, "No exact sketch/mirror pair", ha="center", va="center", transform=ax.transAxes)
    _mark_failure_transitions(ax, report)
    _finish(fig, ax, "Sketch/mirror terminal-current sensitivity", "sketch / mirror ratio")
    entries["topology.png"] = _save_figure(
        fig,
        out_dir / "topology.png",
        series_identities=topology_series,
        failure_markers=failure_markers,
    )
    if set(entries) != set(FIGURE_NAMES):
        raise ValueError("figure renderer did not produce the exact required manifest")
    return {"schema": FIGURE_SCHEMA, "figures": {name: entries[name] for name in FIGURE_NAMES}}

def _write_csv(path: Path, report: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for row in report["checkpoints"]:
        rows.append({"classification": row["classification"], "branch_classification": row["branch_classification"], "branch_threshold_version": row["branch_threshold_version"], "solver": "both", "topology": row["topology"], "bias_V": row["bias_V"], "terminal_current_ratio": row["terminal_current_ratio"]["value"], "maximum_field_ratio": row["maximum_field_ratio"]["value"], "native_source_ratio": row["native_source_ratio"]["value"], "reason": "exact common checkpoint"})
    for row in report["side_only_checkpoints"]:
        rows.append({"classification": "side_only", "branch_classification": None, "branch_threshold_version": report["branch_threshold_version"], "solver": row["solver"], "topology": row["topology"], "bias_V": row["bias_V"], "terminal_current_ratio": None, "maximum_field_ratio": None, "native_source_ratio": None, "reason": "no accepted exact checkpoint from the other solver"})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["classification", "branch_classification", "branch_threshold_version", "solver", "topology", "bias_V", "terminal_current_ratio", "maximum_field_ratio", "native_source_ratio", "reason"])
        writer.writeheader(); writer.writerows(rows)


def _markdown(report: dict[str, Any]) -> str:
    deepest = report["deepest_common_bias_V"]
    deepest_text = str(deepest["value"]) if deepest["classification"] == "available" else f"unavailable ({deepest['reason']})"
    lines = ["# PN2D minimal6 diagnostic sweep comparison", "", DISCLAIMER, "", "## Exact-checkpoint result", "", f"- Deepest common accepted bias: {deepest_text}.", f"- Common exact checkpoints: {len(report['checkpoints'])}.", f"- Recorded rejected transitions: {len(report['failure_transitions'])}.", "- Interpolation is forbidden; solver tails and physical breakdown voltage are not extrapolated.", "", "## Fixed-state recheck", ""]
    lines[10:10] = [f"- {row['topology']} {row['bias_V']:.0f} V: {row['branch_classification']} ({row['branch_threshold_version']})." for row in report["checkpoints"]]
    lines.extend(f"- {row['bias_V']:.0f} V: {row['status']} - {row['reason']}" for row in report["fixed_state_recheck"])
    lines.extend(["", "## Termination", ""])
    lines.extend(f"- {row.get('solver')} {row.get('topology')} {row.get('start_bias_V')} V to {row.get('target_bias_V')} V: {row.get('incomplete_reason', 'rejected transition')}" for row in report["failure_transitions"])
    return "\n".join(lines) + "\n"


def verify_comparison_artifacts(report_path: Path) -> bool:
    """Verify the package and independently rederive its semantic comparison."""
    report_path = Path(report_path)
    if report_path.name != "sweep_comparison.json" or not report_path.is_file():
        raise ValueError("comparison report must be a regular file named sweep_comparison.json")

    def parse_evidence(path: Path, label: str) -> dict[str, Any]:
        def reject_constant(value: str) -> None:
            raise ValueError(f"non-finite JSON constant {value}")
        try:
            parsed = json.loads(
                path.read_text(encoding="utf-8"), parse_constant=reject_constant,
            )
        except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
            raise ValueError(f"{label} is not valid JSON evidence") from exc
        if not isinstance(parsed, dict):
            raise ValueError(f"{label} must be a JSON object")
        return parsed

    report = parse_evidence(report_path, "comparison report")
    validate_bv_comparison_v1(report)
    hashes = report["artifact_hashes"]
    if not isinstance(hashes, dict) or set(hashes) != EXPECTED_ARTIFACT_HASH_NAMES:
        raise ValueError("artifact_hashes must contain the exact package artifact set")
    for name, digest in hashes.items():
        if not _valid_sha256(digest):
            raise ValueError(f"artifact hash is invalid: {name}")
    input_artifacts = report["input_artifacts"]
    if not isinstance(input_artifacts, dict) or set(input_artifacts) != INPUT_ARTIFACT_NAMES:
        raise ValueError("input_artifacts must contain the exact semantic input set")
    input_paths: dict[str, Path] = {}
    inputs: dict[str, dict[str, Any]] = {}
    for name in INPUT_ARTIFACT_NAMES:
        item = input_artifacts[name]
        if not isinstance(item, dict) or set(item) != {"path", "sha256"}:
            raise ValueError(f"input artifact {name} has an invalid contract")
        if not isinstance(item["path"], str) or not _valid_sha256(item["sha256"]):
            raise ValueError(f"input artifact {name} has invalid path/hash fields")
        path = Path(item["path"])
        if not path.is_absolute() or path.resolve() != path or not path.is_file():
            raise ValueError(f"input artifact path must resolve to an absolute regular file: {name}")
        if _sha_bytes(path.read_bytes()) != item["sha256"]:
            raise ValueError(f"input artifact hash mismatch: {name}")
        if hashes[f"input:{name}"] != item["sha256"]:
            raise ValueError(f"input artifact is not hash-addressed: {name}")
        input_paths[name] = path
        inputs[name] = parse_evidence(path, f"input artifact {name}")

    _validate_bound_sweep_input(
        inputs["vela_manifest"], input_paths["vela_manifest"], "vela"
    )
    _validate_bound_sweep_input(
        inputs["sentaurus_manifest"], input_paths["sentaurus_manifest"], "sentaurus"
    )
    validate_formula_difference_v1(inputs["fixed_state_report"])
    expected_report = compare_sweeps(
        inputs["vela_manifest"], inputs["sentaurus_manifest"],
        fixed_state_report=inputs["fixed_state_report"],
    )
    report_semantics = {
        key: value for key, value in report.items() if key not in POST_GENERATION_FIELDS
    }
    expected_semantics = {
        key: value for key, value in expected_report.items()
        if key not in POST_GENERATION_FIELDS
    }
    if report_semantics != expected_semantics:
        raise ValueError(
            "comparison report semantics do not match independently rederived inputs"
        )

    for name in GENERATED_ARTIFACT_NAMES:
        digest = hashes[name]
        path = report_path.parent / name
        if not path.is_file() or _sha_bytes(path.read_bytes()) != digest:
            raise ValueError(f"generated artifact hash mismatch: {name}")
    contract = report.get("figure_contract")
    if not isinstance(contract, dict) or set(contract) != {"schema", "figures"} or contract.get("schema") != FIGURE_SCHEMA:
        raise ValueError("figure contract is invalid")
    entries = contract.get("figures")
    if not isinstance(entries, dict) or set(entries) != set(FIGURE_NAMES):
        raise ValueError("figure contract requires the exact figure manifest")
    expected_entry_keys = {
        "sha256", "width_px", "height_px", "metadata",
        "series_identities", "failure_transition_markers",
    }
    for name in FIGURE_NAMES:
        entry = entries[name]
        if not isinstance(entry, dict) or set(entry) != expected_entry_keys:
            raise ValueError(f"figure contract entry is invalid: {name}")
        path = report_path.parent / name
        if entry["sha256"] != hashes.get(name) or not path.is_file():
            raise ValueError(f"figure hash contract mismatch: {name}")
        payload = path.read_bytes()
        if _sha_bytes(payload) != entry["sha256"] or not payload.startswith(PNG_SIGNATURE):
            raise ValueError(f"figure PNG bytes mismatch: {name}")
        if (
            entry["width_px"] != FIGURE_WIDTH_PX
            or entry["height_px"] != FIGURE_HEIGHT_PX
            or entry["width_px"] < 640
            or entry["height_px"] < 360
            or entry["metadata"] != FIGURE_METADATA
        ):
            raise ValueError(f"figure contract metadata is invalid: {name}")
        try:
            with Image.open(path) as image:
                image.verify()
            with Image.open(path) as image:
                image.load()
                if image.format != "PNG" or image.size != (entry["width_px"], entry["height_px"]):
                    raise ValueError(f"figure PNG dimensions mismatch: {name}")
                if not any(low < high for low, high in image.convert("RGB").getextrema()):
                    raise ValueError(f"figure PNG is blank: {name}")
                for key, value in FIGURE_METADATA.items():
                    if image.info.get(key) != value:
                        raise ValueError(f"figure metadata mismatch: {name}")
                try:
                    series = json.loads(image.info["SeriesIdentities"])
                    markers = json.loads(image.info["FailureTransitionMarkers"])
                except (KeyError, TypeError, json.JSONDecodeError) as exc:
                    raise ValueError(f"figure metadata JSON mismatch: {name}") from exc
                if series != entry["series_identities"] or markers != entry["failure_transition_markers"]:
                    raise ValueError(f"figure metadata identity mismatch: {name}")
        except (OSError, UnidentifiedImageError) as exc:
            raise ValueError(f"figure PNG decode mismatch: {name}") from exc
    with tempfile.TemporaryDirectory(prefix="pn2d-minimal6-figure-verify-") as temp:
        _render_figures(Path(temp), report)
        for name in FIGURE_NAMES:
            actual_path = report_path.parent / name
            expected_path = Path(temp) / name
            with Image.open(actual_path) as actual_image, Image.open(expected_path) as expected_image:
                actual_pixels = actual_image.convert("RGB").tobytes()
                expected_pixels = expected_image.convert("RGB").tobytes()
            if actual_pixels != expected_pixels:
                raise ValueError(
                    f"figure deterministic pixel rerender mismatch: {name}"
                )

    return True


def _preflight_input_artifacts(
    input_artifacts: dict[str, Path] | None,
    expected_objects: dict[str, Any],
) -> dict[str, dict[str, str]]:
    required = set(expected_objects)
    if not isinstance(input_artifacts, dict) or set(input_artifacts) != required:
        actual = set(input_artifacts) if isinstance(input_artifacts, dict) else set()
        raise ValueError(
            f"input_artifacts must have exact semantic keys; expected={sorted(required)}, actual={sorted(actual)}"
        )
    records: dict[str, dict[str, str]] = {}
    for name, path_value in input_artifacts.items():
        if name not in expected_objects:
            raise ValueError(f"input artifact has unknown semantic identity: {name}")
        path = Path(path_value)
        try:
            def reject_constant(value: str) -> None:
                raise ValueError(f"non-finite JSON constant {value}")
            parsed = json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)
        except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
            raise ValueError(f"input artifact {name} is not valid canonical JSON evidence") from exc
        if name == "fixed_state_report":
            if not parsed:
                raise ValueError("fixed_state_report input artifact must be a nonempty schema-valid Phase A report")
            validate_formula_difference_v1(parsed)
        if _canonical_sha(parsed) != _canonical_sha(expected_objects[name]):
            raise ValueError(f"input artifact content does not match supplied object: {name}")
        records[name] = {"path": str(path.resolve()), "sha256": _sha_bytes(path.read_bytes())}
    return records


def write_comparison_package(out_dir: Path, vela_manifest: dict[str, Any], sentaurus_manifest: dict[str, Any], *, fixed_state_report: dict[str, Any], input_artifacts: dict[str, Path] | None = None) -> dict[str, Any]:
    if fixed_state_report:
        validate_formula_difference_v1(fixed_state_report)
    input_records = _preflight_input_artifacts(
        input_artifacts,
        {
            "vela_manifest": vela_manifest,
            "sentaurus_manifest": sentaurus_manifest,
            "fixed_state_report": fixed_state_report,
        },
    )
    assert input_artifacts is not None
    _validate_bound_sweep_input(
        vela_manifest, Path(input_artifacts["vela_manifest"]), "vela")
    _validate_bound_sweep_input(
        sentaurus_manifest, Path(input_artifacts["sentaurus_manifest"]), "sentaurus")
    out_dir.mkdir(parents=True, exist_ok=True)
    report = compare_sweeps(vela_manifest, sentaurus_manifest, fixed_state_report=fixed_state_report)
    csv_path = out_dir / "sweep_comparison.csv"; _write_csv(csv_path, report)
    md_path = out_dir / "sweep_comparison.md"; md_path.write_text(_markdown(report), encoding="utf-8")
    figures = _render_figures(out_dir, report)
    report["figure_contract"] = figures
    hashes = {"sweep_comparison.csv": _sha_bytes(csv_path.read_bytes()), "sweep_comparison.md": _sha_bytes(md_path.read_bytes())}
    hashes.update({name: entry["sha256"] for name, entry in figures["figures"].items()})
    hashes.update({f"input:{name}": item["sha256"] for name, item in input_records.items()})
    report["artifact_hashes"] = hashes
    report["input_artifacts"] = input_records
    validate_bv_comparison_v1(report)
    report_path = out_dir / "sweep_comparison.json"
    report_path.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    verify_comparison_artifacts(report_path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vela-manifest", type=Path, required=True); parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--fixed-state-report", type=Path, required=True); parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    vela = json.loads(args.vela_manifest.read_text(encoding="utf-8")); sentaurus = json.loads(args.sentaurus_manifest.read_text(encoding="utf-8")); fixed = json.loads(args.fixed_state_report.read_text(encoding="utf-8"))
    write_comparison_package(args.out_dir.resolve(), vela, sentaurus, fixed_state_report=fixed, input_artifacts={"vela_manifest": args.vela_manifest, "sentaurus_manifest": args.sentaurus_manifest, "fixed_state_report": args.fixed_state_report})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())