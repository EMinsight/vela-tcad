#!/usr/bin/env python3
"""Contracts shared by the general-Tri3 avalanche diagnostics."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any


SCHEMA_ID = "pn2d_general_tri3_element_edge_avalanche/v1"
EXACT_BIASES_V = (-1.0, -10.0, -20.0)
SENTAURUS_RELEASE = "O-2018.06-SP2"

REQUIRED_RECORD_FIELDS = frozenset(
    {
        "case",
        "topology",
        "bias_V",
        "carrier",
        "vela_cell_id",
        "sentaurus_region_cell_id",
        "node_ids",
        "local_edge",
        "edge_start",
        "edge_end",
        "edge_orientation",
        "contact_adjacent",
        "interior_element",
        "triangle_angles_deg",
        "signed_area_um2",
        "cell_orientation",
        "read_coefficient",
        "read_measure_um2",
        "electric_field_source_sha256",
        "qfp_gradient_source_sha256",
        "low_field_mobility_source_sha256",
        "final_mobility_source_sha256",
        "coefficient_driving_force",
        "current_density_approximation",
        "observation_label",
        "support_status",
        "unit",
        "absolute_error_dex",
    }
)

ALLOWED_OBSERVATION_LABELS = frozenset(
    {
        "native_node",
        "native_element",
        "native_currentplot_integral",
        "box_operator_reconstruction",
        "vela_recomputation",
        "unsupported_native_edge",
    }
)
ALLOWED_SUPPORT_STATUSES = frozenset(
    {"valid", "finite", "below_floor", "geometric_zero", "unsupported"}
)
ALLOWED_DRIVERS = frozenset({"electric_field", "quasi_fermi_gradient"})
ALLOWED_CURRENT_APPROXIMATIONS = frozenset(
    {"element_edge_sg", "aval_dens_grad_qf"}
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_record(record: Mapping[str, Any]) -> None:
    missing = sorted(REQUIRED_RECORD_FIELDS - set(record))
    _require(not missing, f"missing required record fields: {missing}")
    _require(record["carrier"] in {"electron", "hole"}, "invalid carrier")
    _require(
        record["observation_label"] in ALLOWED_OBSERVATION_LABELS,
        "invalid observation label",
    )
    _require(
        record["support_status"] in ALLOWED_SUPPORT_STATUSES,
        "invalid support status",
    )
    _require(
        record["coefficient_driving_force"] in ALLOWED_DRIVERS,
        "invalid coefficient driving force",
    )
    _require(
        record["current_density_approximation"]
        in ALLOWED_CURRENT_APPROXIMATIONS,
        "invalid current-density approximation",
    )
    _require(
        record["edge_orientation"] in {-1, 1},
        "invalid carrier or edge orientation",
    )
    _require(
        bool(record["contact_adjacent"]) != bool(record["interior_element"]),
        "contact/interior classification mismatch",
    )
    _require(
        len(record["triangle_angles_deg"]) == 3,
        "triangle angle count mismatch",
    )
    _require(
        abs(sum(float(value) for value in record["triangle_angles_deg"]) - 180.0)
        <= 1.0e-8,
        "triangle angle sum mismatch",
    )
    _require(
        record["cell_orientation"] in {"ccw", "cw"},
        "invalid cell orientation",
    )
    _require(
        float(record["signed_area_um2"]) != 0.0,
        "degenerate triangle",
    )

    if record["observation_label"] == "native_element":
        _require(
            record["current_density_approximation"] != "element_edge_sg",
            "reconstructed edge current mislabeled as native",
        )
    if record["support_status"] == "geometric_zero":
        _require(
            float(record["read_coefficient"]) == 0.0,
            "geometric zero has nonzero coefficient",
        )
        _require(
            record["absolute_error_dex"] in {"", None},
            "geometric zero converted to finite dex",
        )


def validate_cell_mapping(
    records: Iterable[Mapping[str, Any]],
    expected: Mapping[tuple[str, int], int],
) -> None:
    for record in records:
        key = (str(record["topology"]), int(record["vela_cell_id"]))
        _require(key in expected, f"missing expected cell mapping for {key}")
        _require(
            int(record["sentaurus_region_cell_id"]) == int(expected[key]),
            f"wrong cell permutation for {key}",
        )


def validate_driver_contract(record: Mapping[str, Any]) -> None:
    global_default = record.get("global_default_driver")
    effective = record.get("effective_driver")
    _require(
        global_default == "quasi_fermi_gradient",
        "global Sentaurus default must remain quasi_fermi_gradient",
    )
    if bool(record["contact_adjacent"]) and not bool(
        record.get("use_quasi_fermi_at_contacts", False)
    ):
        _require(
            effective == "electric_field",
            "contact fallback mislabeled",
        )
    elif record["coefficient_driving_force"] == "quasi_fermi_gradient":
        _require(
            effective == "quasi_fermi_gradient",
            "interior or forced-QFP driver mismatch",
        )


def validate_source_manifests(
    manifests: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, str]]:
    _require(bool(manifests), "no source manifests")
    shared_hashes: dict[str, dict[str, str]] = {}
    for label, manifest in manifests.items():
        _require(manifest.get("schema") == SCHEMA_ID, f"{label}: schema mismatch")
        _require(
            manifest.get("sentaurus_release") == SENTAURUS_RELEASE,
            "Sentaurus releases differ",
        )
        _require(manifest.get("status") == "passed", f"{label}: status mismatch")
        biases = tuple(float(value) for value in manifest.get("exact_biases_V", ()))
        _require(biases == EXACT_BIASES_V, f"{label}: bias matrix mismatch")
        case_hashes = manifest.get("case_hashes")
        _require(isinstance(case_hashes, Mapping), f"{label}: missing case hashes")
        for case, raw_hashes in case_hashes.items():
            _require(isinstance(raw_hashes, Mapping), f"{label}: invalid case hashes")
            hashes = {str(key): str(value) for key, value in raw_hashes.items()}
            _require(
                {"tdr", "models.par"} <= set(hashes),
                f"{label}: incomplete static bundle hashes",
            )
            previous = shared_hashes.setdefault(str(case), hashes)
            _require(
                previous == hashes,
                f"{label}: static bundle hash mismatch for {case}",
            )
    return shared_hashes
