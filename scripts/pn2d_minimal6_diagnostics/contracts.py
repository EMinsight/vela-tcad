"""Typed identities shared by PN2D Minimal6 diagnostic analyses."""

from dataclasses import dataclass
from enum import StrEnum
import math
import re


class SourceKind(StrEnum):
    SENTAURUS = "sentaurus"
    VELA = "vela"
    DERIVED = "derived"


class SupportKind(StrEnum):
    NODE = "node"
    EDGE = "edge"
    CELL = "cell"
    TERMINAL = "terminal"


class BranchKind(StrEnum):
    LEAKAGE_LIKE = "leakage_like"
    MULTIPLICATION_LIKE = "multiplication_like"
    UNIDENTIFIED = "unidentified"


_SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")


@dataclass(frozen=True)
class StateIdentity:
    run_id: str
    topology: str
    bias_V: float

    def __post_init__(self) -> None:
        if not isinstance(self.run_id, str) or not self.run_id:
            raise ValueError("state identity requires a run_id")
        if not isinstance(self.topology, str) or not self.topology:
            raise ValueError("state identity requires a topology")
        if isinstance(self.bias_V, bool) or not isinstance(self.bias_V, (int, float)):
            raise TypeError("state bias must be numeric")
        if not math.isfinite(float(self.bias_V)):
            raise ValueError("state bias must be finite")


@dataclass(frozen=True)
class QuantityRecord:
    state: StateIdentity
    carrier: str
    support_kind: SupportKind
    support_id: str
    quantity: str
    source: SourceKind
    formula_version: str
    value: float
    unit: str
    sign_convention: str
    raw_source_path: str
    raw_source_sha256: str
    geometric_zero_reason: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.state, StateIdentity):
            raise TypeError("record state must be a StateIdentity")
        if not isinstance(self.support_kind, SupportKind):
            raise TypeError("record support_kind must be a SupportKind")
        if not isinstance(self.source, SourceKind):
            raise TypeError("record source must be a SourceKind")
        required = {
            "carrier": self.carrier,
            "support_id": self.support_id,
            "quantity": self.quantity,
            "formula_version": self.formula_version,
            "unit": self.unit,
            "sign_convention": self.sign_convention,
            "raw_source_path": self.raw_source_path,
        }
        if any(not isinstance(value, str) or not value for value in required.values()):
            raise ValueError(f"record requires non-empty {', '.join(required)}")
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise TypeError("record value must be numeric")
        if not math.isfinite(float(self.value)):
            raise ValueError("record value must be finite")
        if not isinstance(self.raw_source_sha256, str) or not _SHA256.fullmatch(
            self.raw_source_sha256
        ):
            raise ValueError("record raw source hash must be a SHA-256 hex digest")
        if self.geometric_zero_reason is not None and (
            not isinstance(self.geometric_zero_reason, str) or not self.geometric_zero_reason
        ):
            raise ValueError("geometric-zero reason must be a non-empty string")

    @property
    def key(self) -> tuple[object, ...]:
        return (
            self.state.run_id,
            self.state.topology,
            float(self.state.bias_V),
            self.carrier,
            self.support_kind,
            self.support_id,
            self.quantity,
            self.source,
            self.formula_version,
        )


def classify_pair(
    reference_current_A_per_um: float,
    candidate_current_A_per_um: float,
    multiplication_floor_A_per_um: float = 1.0e-9,
) -> BranchKind:
    values = (reference_current_A_per_um, candidate_current_A_per_um)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in values):
        return BranchKind.UNIDENTIFIED
    if not all(math.isfinite(float(value)) for value in values):
        return BranchKind.UNIDENTIFIED
    if not math.isfinite(multiplication_floor_A_per_um) or multiplication_floor_A_per_um <= 0.0:
        raise ValueError("multiplication floor must be finite and positive")
    if reference_current_A_per_um == 0.0 or candidate_current_A_per_um == 0.0:
        return BranchKind.UNIDENTIFIED
    if max(abs(reference_current_A_per_um), abs(candidate_current_A_per_um)) >= multiplication_floor_A_per_um:
        return BranchKind.MULTIPLICATION_LIKE
    return BranchKind.LEAKAGE_LIKE
