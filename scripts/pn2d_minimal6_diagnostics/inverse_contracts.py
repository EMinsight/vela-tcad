from dataclasses import dataclass
from enum import Enum
import math
from typing import Any


class SupportKind(str, Enum):
    NODE = "node"
    EDGE = "edge"
    CELL = "cell"
    CONTACT = "contact"
    INTEGRATED = "integrated"


class SampleStatus(str, Enum):
    VALID = "valid"
    GEOMETRIC_ZERO = "geometric_zero"
    BELOW_FLOOR = "below_numerical_floor"
    MISSING_FIELD = "missing_field"
    INCOMPATIBLE_SUPPORT = "incompatible_support"
    INVALID_UNIT = "invalid_unit"
    DIRECTION_UNDEFINED = "direction_undefined"
    BRANCH_AMBIGUOUS = "coefficient_branch_ambiguous"
    EXPONENTIAL_UNDERFLOW = "exponential_underflow"
    NONFINITE = "nonfinite"


class Identifiability(str, Enum):
    IDENTIFIED = "identified"
    CONSISTENT_NONUNIQUE = "consistent_nonunique"
    CONFOUNDED = "confounded"
    INSUFFICIENT_DATA = "insufficient_data"
    REJECTED = "rejected"


@dataclass(frozen=True)
class Observation:
    solver: str
    topology: str
    bias_V: float
    support_kind: SupportKind
    support_id: int | str
    quantity: str
    component: str
    raw_value: float | None
    raw_unit: str
    value_si: float | None
    unit_si: str
    coordinate_frame: str
    orientation: str
    conversion: str
    status: SampleStatus
    source_path: str
    source_sha256: str

    @property
    def key(self) -> tuple[str, str, float, str, int | str, str, str]:
        return (
            self.solver,
            self.topology,
            self.bias_V,
            self.support_kind.value,
            self.support_id,
            self.quantity,
            self.component,
        )


@dataclass(frozen=True)
class AcceptanceThresholds:
    field_median_relative: float = 0.02
    field_median_angle_deg: float = 1.0
    gradient_median_abs_dex: float = 0.1
    gradient_p95_abs_dex: float = 0.3
    gradient_median_angle_deg: float = 5.0
    integrated_generation_abs_dex: float = 0.1
    local_generation_abs_dex: float = 0.3
    replacement_closure_abs_dex: float = 1.0e-10


@dataclass(frozen=True)
class CandidateMetric:
    candidate: str
    quantity: str
    carrier: str
    split: str
    topology: str
    bias_V: float | None
    support_kind: SupportKind
    valid_count: int
    median_abs_error: float | None
    p95_abs_error: float | None
    median_angle_deg: float | None
    classification: Identifiability


_REPORT_KEYS = {"schema", "diagnostic_only", "phase_base", "payload"}
_PAYLOAD_KEYS = {
    "input_manifest_sha256",
    "discovery_keys",
    "holdout_keys",
    "thresholds",
    "field_inventory",
    "sample_status_counts",
    "candidate_metrics",
    "classifications",
    "replacement_closure",
    "localization_control",
    "sentaurus_version",
    "production_cpp_changed",
}


def validate_inverse_report_v1(report: dict[str, Any]) -> dict[str, Any]:
    if set(report) != _REPORT_KEYS:
        raise ValueError("inverse report top-level contract mismatch")
    if report["schema"] != "vela.pn2d_minimal6_physics_inverse_audit.v1":
        raise ValueError("inverse report schema mismatch")
    if report["diagnostic_only"] is not True or report["phase_base"] != "a5524cf":
        raise ValueError("inverse report provenance mismatch")

    payload = report["payload"]
    if not isinstance(payload, dict) or set(payload) != _PAYLOAD_KEYS:
        raise ValueError("inverse report payload contract mismatch")

    allowed = {item.value for item in Identifiability}
    for row in payload["classifications"]:
        if row.get("classification") not in allowed:
            raise ValueError("unknown inverse classification")
    for row in payload["candidate_metrics"]:
        if row.get("classification") not in allowed:
            raise ValueError("unknown inverse classification")
        for name in ("median_abs_error", "p95_abs_error", "median_angle_deg"):
            value = row.get(name)
            if value is not None and not math.isfinite(float(value)):
                raise ValueError("non-finite inverse metric")
    return report


def classify_numeric_sample(
    value: float | None, *, floor: float, geometric_zero: bool = False
) -> SampleStatus:
    if value is None:
        return SampleStatus.MISSING_FIELD
    if not math.isfinite(value):
        return SampleStatus.NONFINITE
    if geometric_zero and value == 0.0:
        return SampleStatus.GEOMETRIC_ZERO
    if abs(value) < floor:
        return SampleStatus.BELOW_FLOOR
    return SampleStatus.VALID
