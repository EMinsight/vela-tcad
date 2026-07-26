#!/usr/bin/env python3
"""Contracts for PN2D high-bias process and source-Jacobian evidence."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from typing import Any


SCHEMA_ID = "vela.pn2d_high_bias_process_jacobian.v1"
SENTAURUS_RELEASE = "O-2018.06-SP2"
EXACT_HIGH_BIAS_V = (
    -10.0,
    -18.0,
    -19.0,
    -19.5,
    -19.8,
    -19.9,
    -19.95,
    -20.0,
)
REQUIRED_DERIVATIVE_STAGES = (
    "carrier_statistics",
    "sg_bernoulli",
    "low_field_mobility",
    "high_field_mobility",
    "current_vector_reconstruction",
    "current_abs_branch",
    "avalanche_alpha",
    "element_source_distribution",
    "node_accumulation",
    "continuity_scaling",
)

REQUIRED_RECORD_FIELDS = frozenset(
    {
        "schema",
        "topology",
        "bias_V",
        "carrier",
        "residual_variable",
        "state_variable",
        "process_stage",
        "analytic_contribution",
        "fd_contribution",
        "analytic_total",
        "fd_total",
        "derivative_status",
        "observation_label",
        "observation_provenance",
        "support_status",
        "error_dex",
        "unit",
        "fd_steps_V",
        "residual_config_sha256",
        "jacobian_config_sha256",
        "state_sha256",
        "deck_sha256",
        "tdr_sha256",
        "mesh_sha256",
        "parameters_sha256",
        "sentaurus_release",
    }
)
REQUIRED_INPUT_HASHES = frozenset({"deck", "tdr", "mesh", "parameters", "state"})
ALLOWED_CARRIERS = frozenset({"electron", "hole"})
ALLOWED_RESIDUAL_VARIABLES = frozenset(
    {"poisson", "electron_continuity", "hole_continuity"}
)
ALLOWED_STATE_VARIABLES = frozenset({"psi", "electron_qfp", "hole_qfp"})
ALLOWED_DERIVATIVE_STATUSES = frozenset(
    {"finite", "zero", "below_floor", "nonsmooth", "unsupported", "invalid"}
)
ALLOWED_SUPPORT_STATUSES = frozenset(
    {"valid", "finite", "zero", "below_floor", "nonsmooth", "unsupported", "invalid"}
)
ALLOWED_OBSERVATION_LABELS = frozenset(
    {
        "native_node",
        "native_element",
        "native_currentplot_integral",
        "native_directed_edge",
        "operator_replay",
        "box_operator_reconstruction",
        "unsupported_native_edge",
        "unsupported",
    }
)
ALLOWED_OBSERVATION_PROVENANCE = frozenset(
    {"native", "reconstructed", "unsupported"}
)
HASH_RE = re.compile(r"^[0-9a-f]{64}$")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_hash(value: Any, label: str) -> None:
    _require(
        isinstance(value, str) and HASH_RE.fullmatch(value) is not None,
        f"invalid {label} SHA-256",
    )


def true_relative_error(analytic: float, finite_difference: float) -> float:
    """Return true relative error without an absolute floor of one."""

    denominator = max(abs(float(analytic)), abs(float(finite_difference)))
    if denominator == 0.0:
        return 0.0
    return abs(float(analytic) - float(finite_difference)) / denominator


def validate_process_record(record: Mapping[str, Any]) -> None:
    missing = sorted(REQUIRED_RECORD_FIELDS - set(record))
    _require(not missing, f"missing process record fields: {missing}")
    _require(record["schema"] == SCHEMA_ID, "process record schema mismatch")
    _require(
        float(record["bias_V"]) in EXACT_HIGH_BIAS_V,
        "non-exact bias in process record",
    )
    _require(record["carrier"] in ALLOWED_CARRIERS, "invalid carrier")
    _require(
        record["residual_variable"] in ALLOWED_RESIDUAL_VARIABLES,
        "invalid residual variable",
    )
    _require(
        record["state_variable"] in ALLOWED_STATE_VARIABLES,
        "invalid state variable",
    )
    _require(
        record["process_stage"] in REQUIRED_DERIVATIVE_STAGES,
        "invalid derivative process stage",
    )
    _require(
        record["derivative_status"] in ALLOWED_DERIVATIVE_STATUSES,
        "invalid derivative status",
    )
    _require(
        record["support_status"] in ALLOWED_SUPPORT_STATUSES,
        "invalid support status",
    )
    provenance = record["observation_provenance"]
    label = record["observation_label"]
    _require(
        provenance in ALLOWED_OBSERVATION_PROVENANCE,
        "invalid observation provenance",
    )
    if provenance == "reconstructed" and str(label).startswith("native"):
        raise ValueError("reconstructed current mislabeled as native")
    _require(label in ALLOWED_OBSERVATION_LABELS, "invalid observation label")
    if label in {"operator_replay", "box_operator_reconstruction"}:
        _require(
            provenance == "reconstructed",
            "operator replay must be labeled reconstructed",
        )
    if label == "unsupported_native_edge":
        _require(
            provenance == "unsupported",
            "unsupported native edge observation provenance mismatch",
        )
    if record["support_status"] == "zero":
        _require(
            record["error_dex"] in {None, ""},
            "zero converted to finite dex",
        )

    for name in (
        "analytic_contribution",
        "fd_contribution",
        "analytic_total",
        "fd_total",
    ):
        _require(math.isfinite(float(record[name])), f"non-finite {name}")
    steps = tuple(float(value) for value in record["fd_steps_V"])
    _require(len(steps) == 3, "finite-difference step count mismatch")
    _require(all(math.isfinite(step) and step > 0.0 for step in steps), "invalid FD step")
    _require(
        steps[0] > steps[1] > steps[2],
        "finite-difference steps must be strictly descending",
    )

    for name in (
        "residual_config_sha256",
        "jacobian_config_sha256",
        "state_sha256",
        "deck_sha256",
        "tdr_sha256",
        "mesh_sha256",
        "parameters_sha256",
    ):
        _require_hash(record[name], name)
    _require(
        record["residual_config_sha256"] == record["jacobian_config_sha256"],
        "residual/Jacobian configuration mismatch",
    )
    _require(
        record["sentaurus_release"] == SENTAURUS_RELEASE,
        "Sentaurus release mismatch",
    )


def validate_derivative_lattice(records: Iterable[Mapping[str, Any]]) -> None:
    groups: dict[tuple[Any, ...], list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        validate_process_record(record)
        key = (
            record["topology"],
            float(record["bias_V"]),
            record["carrier"],
            record["residual_variable"],
            record["state_variable"],
        )
        groups[key].append(record)
    _require(bool(groups), "empty derivative lattice")

    expected = set(REQUIRED_DERIVATIVE_STAGES)
    for key, group in groups.items():
        stages = [str(record["process_stage"]) for record in group]
        missing = sorted(expected - set(stages))
        _require(
            not missing,
            f"missing derivative contribution for {key}: {missing}",
        )
        _require(
            len(stages) == len(set(stages)),
            f"duplicate derivative contribution for {key}",
        )
        _require(
            set(stages) == expected,
            f"unexpected derivative contribution for {key}",
        )
        analytic_total_values = {float(record["analytic_total"]) for record in group}
        fd_total_values = {float(record["fd_total"]) for record in group}
        _require(
            len(analytic_total_values) == 1 and len(fd_total_values) == 1,
            f"inconsistent derivative totals for {key}",
        )
        analytic_total = next(iter(analytic_total_values))
        fd_total = next(iter(fd_total_values))
        analytic_sum = math.fsum(
            float(record["analytic_contribution"]) for record in group
        )
        fd_sum = math.fsum(float(record["fd_contribution"]) for record in group)
        _require(
            true_relative_error(analytic_sum, analytic_total) <= 1.0e-12,
            f"analytic derivative contribution closure failed for {key}",
        )
        _require(
            true_relative_error(fd_sum, fd_total) <= 1.0e-12,
            f"FD derivative contribution closure failed for {key}",
        )


def _validate_manifest(manifest: Mapping[str, Any]) -> None:
    _require(manifest.get("schema") == SCHEMA_ID, "manifest schema mismatch")
    _require(
        manifest.get("status") == "red_contract_frozen",
        "manifest status mismatch",
    )
    _require(
        manifest.get("sentaurus_release") == SENTAURUS_RELEASE,
        "Sentaurus release mismatch",
    )
    biases = tuple(float(value) for value in manifest.get("exact_biases_V", ()))
    _require(biases == EXACT_HIGH_BIAS_V, "exact bias lattice mismatch")
    input_hashes = manifest.get("input_hashes")
    _require(isinstance(input_hashes, Mapping), "missing input hashes")
    _require(
        set(input_hashes) == REQUIRED_INPUT_HASHES,
        "input hash key set mismatch",
    )
    for name, value in input_hashes.items():
        _require_hash(value, str(name))
    _require_hash(manifest.get("residual_config_sha256"), "residual configuration")
    _require_hash(manifest.get("jacobian_config_sha256"), "Jacobian configuration")
    _require(
        manifest["residual_config_sha256"] == manifest["jacobian_config_sha256"],
        "residual/Jacobian configuration mismatch",
    )


def validate_manifest_pair(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
) -> None:
    _validate_manifest(first)
    _validate_manifest(second)
    _require(
        dict(first["input_hashes"]) == dict(second["input_hashes"]),
        "input provenance mismatch",
    )
    _require(
        first["residual_config_sha256"] == second["residual_config_sha256"],
        "paired residual configuration mismatch",
    )
    _require(
        first["jacobian_config_sha256"] == second["jacobian_config_sha256"],
        "paired Jacobian configuration mismatch",
    )
