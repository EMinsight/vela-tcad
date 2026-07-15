from dataclasses import dataclass
from enum import StrEnum
import math

class SourceKind(StrEnum): SENTRAURUS="sentaurus"; VELA="vela"; DERIVED="derived"
class SupportKind(StrEnum): NODE="node"; EDGE="edge"; CELL="cell"; TERMINAL="terminal"
class BranchKind(StrEnum): LEAKAGE_LIKE="leakage_like"; MULTIPLICATION_LIKE="multiplication_like"; UNIDENTIFIED="unidentified"
@dataclass(frozen=True)
class StateIdentity: run_id:str; topology:str; bias_V:float
@dataclass(frozen=True)
class QuantityRecord:
    state:StateIdentity; carrier:str; support_kind:SupportKind; support_id:str; quantity:str; source:SourceKind; formula_version:str; value:float; unit:str; sign_convention:str; raw_source_path:str; raw_source_sha256:str; geometric_zero_reason:str|None=None
    def __post_init__(self):
        if not self.state.run_id or not self.quantity or not self.formula_version or not self.sign_convention or not self.raw_source_path or len(self.raw_source_sha256) != 64 or not math.isfinite(self.value): raise ValueError("record requires identifiers, provenance, sign convention, and finite value")
def classify_pair(reference_current_A_per_um:float,candidate_current_A_per_um:float,multiplication_floor_A_per_um:float=1e-9)->BranchKind:
    if not all(math.isfinite(v) for v in (reference_current_A_per_um,candidate_current_A_per_um)): return BranchKind.UNIDENTIFIED
    return BranchKind.MULTIPLICATION_LIKE if max(abs(reference_current_A_per_um),abs(candidate_current_A_per_um))>=multiplication_floor_A_per_um else BranchKind.LEAKAGE_LIKE
