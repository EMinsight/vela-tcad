# PN2D Minimal6 Physics Inverse Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and execute a provenance-complete, report-only inverse audit of PN2D Minimal6 potential, electric field, quasi-Fermi gradients, current densities, avalanche driving force, alpha, and generation without changing production C++ formulas.

**Architecture:** Extend the existing `scripts/pn2d_minimal6_diagnostics` package with focused contract, input, field, transport, avalanche, replacement, and report modules. Consume the sealed Task 8 Vela and Sentaurus sweeps plus a hash-bound supplemental Sentaurus field matrix, normalize them onto explicit node/edge/cell supports, evaluate candidate formulas, and replay one-factor/staged replacements through a declared dependency graph. Emit deterministic machine-readable tables, static figures, a scientific Markdown report, and an independent verifier.

**Tech Stack:** Python 3.11, standard library (`argparse`, `csv`, `dataclasses`, `enum`, `hashlib`, `json`, `math`, `pathlib`, `unittest`), NumPy, Matplotlib, existing Minimal6 topology/diagnostic helpers, existing Sentaurus importer and remote state exporter.

## Global Constraints

- This phase is diagnostic only; do not modify files under `include/` or `src/`.
- Treat `a5524cf` as the production-code baseline and prove that `include/` and `src/` are unchanged from it.
- Use immutable, hash-addressed inputs; never overwrite raw Vela or Sentaurus values.
- Use exact common checkpoints `-1..-20 V` for sketch and mirror; do not interpolate bias.
- Use discovery checkpoints `sketch x {-1,-4,-8,-12,-16,-19,-20 V}` only for ranking or any global constant selection.
- Use mirror plus all remaining exact checkpoints as holdout evidence; never fit by bias, node, edge, cell, carrier, or topology.
- Preserve typed `identified`, `consistent_nonunique`, `confounded`, `insufficient_data`, and `rejected` conclusions.
- Exclude typed geometric zeros, below-floor samples, undefined directions, and non-finite samples from ordinary relative-error statistics; never convert missing data to zero.
- Require electric-field median magnitude error <= 2% and median direction error <= 1 degree for identification.
- Require quasi-Fermi-gradient/current median absolute log10 error <= 0.1 dex, 95th percentile <= 0.3 dex, and median direction error <= 5 degrees.
- Require valid high-field integrated alpha/current generation error <= 0.1 dex and local support error <= 0.3 dex.
- Require finite nonzero replacement closure <= `1e-10 dex`.
- Use the bundled Python executable `C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`; add no dependencies.
- Preserve the existing untracked `docs/validation/figures/pn2d_minimal6_{sketch,mirror}_mesh_nodes.png` files.

## File Structure

- Create `scripts/pn2d_minimal6_diagnostics/inverse_contracts.py`: enums, immutable record types, thresholds, schema validation, and typed sample classification.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_inputs.py`: strict sweep/supplement validation, hash binding, field inventory, and canonical observation loading.
- Modify `scripts/export_pn2d_minimal6_states.py`: support an exact caller-declared bias matrix while retaining the legacy six-state default and record remote Sentaurus version provenance.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_fields.py`: topology transforms, support reconstruction, potential gradients, electric-field candidates, and vector metrics.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_transport.py`: quasi-Fermi gradients, current projections, current-inverted effective gradients, and transport metrics.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_avalanche.py`: Van Overstraeten inversion, driving-force candidates, generation reconstruction, and support integration.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_replacements.py`: dependency graph, one-factor/staged/reverse replacement, closure, ranking, and identifiability classification.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_report.py`: deterministic CSV/JSON/Markdown/manifest writers.
- Create `scripts/pn2d_minimal6_diagnostics/inverse_plots.py`: deterministic PNG/PDF figures.
- Create `scripts/diagnose_pn2d_minimal6_physics_inverse_audit.py`: report CLI.
- Create `scripts/verify_pn2d_minimal6_physics_inverse_audit.py`: independent semantic and hash verifier.
- Create seven focused regression modules under `tests/regression/` matching the tasks below; Task 3 extends the existing state-export regression module.
- Modify `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`: append the authoritative inverse-audit result only after the report passes.

---

### Task 1: Freeze inverse-audit contracts and typed outcomes

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_contracts.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_contracts.py`

**Interfaces:**
- Produces: `SupportKind`, `SampleStatus`, `Identifiability`, `Observation`, `CandidateMetric`, `AcceptanceThresholds`, `classify_numeric_sample()`, and `validate_inverse_report_v1()`.
- Consumes: no new interfaces.

- [ ] **Step 1: Write the failing enum, record, and schema tests**

```python
import math
import unittest

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds, Identifiability, Observation, SampleStatus,
    SupportKind, classify_numeric_sample, validate_inverse_report_v1,
)


class InverseContractsTest(unittest.TestCase):
    def test_numeric_statuses_do_not_turn_missing_into_zero(self):
        self.assertEqual(classify_numeric_sample(None, floor=1e-30), SampleStatus.MISSING_FIELD)
        self.assertEqual(classify_numeric_sample(0.0, floor=1e-30, geometric_zero=True), SampleStatus.GEOMETRIC_ZERO)
        self.assertEqual(classify_numeric_sample(1e-40, floor=1e-30), SampleStatus.BELOW_FLOOR)
        self.assertEqual(classify_numeric_sample(math.inf, floor=1e-30), SampleStatus.NONFINITE)
        self.assertEqual(classify_numeric_sample(2.0, floor=1e-30), SampleStatus.VALID)

    def test_observation_key_is_complete_and_immutable(self):
        row = Observation(
            "sentaurus", "sketch", -12.0, SupportKind.NODE, 1,
            "electric_field", "x", -64036.5, "V*cm^-1", -6.40365e6, "V/m",
            "sentaurus_xy", "global_vector", "multiply_by_100",
            SampleStatus.VALID, "field.csv", "0" * 64,
        )
        self.assertEqual(row.key, ("sentaurus", "sketch", -12.0, "node", 1,
                                   "electric_field", "x"))
        with self.assertRaises(Exception):
            row.value_si = 0.0

    def test_thresholds_and_report_schema_are_exact(self):
        limits = AcceptanceThresholds()
        self.assertEqual(limits.gradient_median_abs_dex, 0.1)
        self.assertEqual(limits.gradient_p95_abs_dex, 0.3)
        report = {
            "schema": "vela.pn2d_minimal6_physics_inverse_audit.v1",
            "diagnostic_only": True, "phase_base": "a5524cf",
            "payload": {
                "input_manifest_sha256": "0" * 64,
                "discovery_keys": [["sketch", -1.0]],
                "holdout_keys": [["mirror", -1.0]], "thresholds": limits.__dict__,
                "field_inventory": {}, "sample_status_counts": {"valid": 1},
                "candidate_metrics": [],
                "classifications": [{"candidate": "triangle_gradient",
                                      "classification": Identifiability.IDENTIFIED.value}],
                "replacement_closure": [], "localization_control": {},
                "sentaurus_version": "O-2018.06-SP2", "production_cpp_changed": False,
            },
        }
        validate_inverse_report_v1(report)
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```powershell
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" -m unittest discover -s tests\regression -p "test_pn2d_minimal6_inverse_contracts.py" -v
```

Expected: FAIL with `ModuleNotFoundError: scripts.pn2d_minimal6_diagnostics.inverse_contracts`.

- [ ] **Step 3: Implement the immutable contracts**

```python
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
        return (self.solver, self.topology, self.bias_V, self.support_kind.value,
                self.support_id, self.quantity, self.component)


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
    "input_manifest_sha256", "discovery_keys", "holdout_keys", "thresholds",
    "field_inventory", "sample_status_counts", "candidate_metrics",
    "classifications", "replacement_closure", "localization_control",
    "sentaurus_version", "production_cpp_changed",
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
        for name in ("median_abs_error", "p95_abs_error", "median_angle_deg"):
            value = row.get(name)
            if value is not None and not math.isfinite(float(value)):
                raise ValueError("non-finite inverse metric")
    return report

def classify_numeric_sample(value: float | None, *, floor: float,
                            geometric_zero: bool = False) -> SampleStatus:
    if value is None:
        return SampleStatus.MISSING_FIELD
    if not math.isfinite(value):
        return SampleStatus.NONFINITE
    if geometric_zero and value == 0.0:
        return SampleStatus.GEOMETRIC_ZERO
    if abs(value) < floor:
        return SampleStatus.BELOW_FLOOR
    return SampleStatus.VALID
```

Define `CandidateMetric` as a frozen dataclass containing candidate, quantity,
carrier, split, topology, bias, support, valid count, median error, p95 error,
median angle, and classification. Make `validate_inverse_report_v1()` fail
closed on extra/missing top-level keys, non-diagnostic reports, a phase base
other than `a5524cf`, unknown classifications, or non-finite metrics.

- [ ] **Step 4: Run the contract tests and the existing Minimal6 schema tests**

Run:

```powershell
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" -m unittest discover -s tests\regression -p "test_pn2d_minimal6_inverse_contracts.py" -v
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" -m unittest discover -s tests\regression -p "test_pn2d_minimal6_diagnostic_contracts.py" -v
```

Expected: both commands PASS.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_contracts.py tests/regression/test_pn2d_minimal6_inverse_contracts.py
D:\msys64\usr\bin\git.exe commit -m "Add Minimal6 inverse audit contracts"
```

### Task 2: Validate and canonicalize the sealed input roots

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_inputs.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_inputs.py`

**Interfaces:**
- Consumes: Task 1 `Observation`, `SampleStatus`, and `SupportKind`.
- Produces: frozen `InputBundle`, `load_input_bundle(vela_root, sentaurus_root, supplemental_root)`, `field_inventory(bundle)`, `canonical_observations(bundle)`, and `write_input_manifest(bundle, path)`.

- [ ] **Step 1: Write failing fixture tests for exact matrices, hashes, fields, and supplement binding**

Create a temporary fixture with two topologies, exact `-1..-20 V` Vela and
Sentaurus records, and a supplemental state manifest carrying the same forty
keys. Assert:

```python
bundle = load_input_bundle(vela_root, sentaurus_root, supplemental_root)
self.assertEqual(len(bundle.common_keys), 40)
self.assertEqual(bundle.discovery_keys,
                 tuple(("sketch", float(v)) for v in (-1, -4, -8, -12, -16, -19, -20)))
self.assertEqual(len(bundle.holdout_keys), 33)
self.assertEqual(field_inventory(bundle)["eMobility"]["unit"], "cm^2*V^-1*s^-1")

(vela_root / "vela/sketch/states/segment_00_bias_m1p000000.csv").write_text("tampered")
with self.assertRaisesRegex(ValueError, "state hash mismatch"):
    load_input_bundle(vela_root, sentaurus_root, supplemental_root)
```

Also test rejection of duplicate/inexact biases, absolute paths escaping a
root, a supplemental topology/bias mismatch, a field component/unit mismatch,
and absent `eMobility` classified as missing rather than zero.

- [ ] **Step 2: Run the input test and verify RED**

Run the `test_pn2d_minimal6_inverse_inputs.py` module with unittest discovery.

Expected: FAIL because `inverse_inputs` does not exist.

- [ ] **Step 3: Implement strict loading and canonical row creation**

Use these exact constants and types:

```python
COMMON_BIASES = tuple(float(-value) for value in range(1, 21))
TOPOLOGIES = ("sketch", "mirror")
DISCOVERY_KEYS = tuple(("sketch", value) for value in (-1.0, -4.0, -8.0,
                                                        -12.0, -16.0, -19.0, -20.0))
REQUIRED_SENTARUS_FIELDS = {
    "ElectrostaticPotential": (1, "V"),
    "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"),
    "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"),
    "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"),
    "hCurrentDensity": (2, "A*cm^-2"),
    "eAlphaAvalanche": (1, "cm^-1"),
    "hAlphaAvalanche": (1, "cm^-1"),
    "ImpactIonization": (1, "cm^-3*s^-1"),
}
SUPPLEMENTAL_FIELDS = {
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eVelocity": (1, "cm*s^-1"),
    "hVelocity": (1, "cm*s^-1"),
}
```

Resolve every declared member beneath its bound root, recompute SHA-256 before
reading, require exact actual/requested bias equality within `1e-12 V`, and
join supplemental fields only on exact `(topology, bias, canonical_node_id)`.
Convert coordinates to metres, density to `m^-3`, field to `V/m`, current to
`A/m^2`, mobility to `m^2/(V*s)`, alpha to `m^-1`, velocity to `m/s`, and keep
generation in `m^-3*s^-1`. Store raw and SI values, both units, coordinate
frame, orientation, and the named conversion in every `Observation`. Copy the
executable hashes and tracked-source hashes bound by the sweep manifests into
the input manifest.

- [ ] **Step 4: Run focused and adversarial input tests**

Run the new input module plus `test_pn2d_minimal6_diagnostic_sweep.py` and
`test_pn2d_minimal6_state_export.py`.

Expected: PASS; every adversarial mutation fails before an output directory is created.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_inputs.py tests/regression/test_pn2d_minimal6_inverse_inputs.py
D:\msys64\usr\bin\git.exe commit -m "Validate Minimal6 inverse audit inputs"
```

### Task 3: Generalize the supplemental Sentaurus field matrix safely

**Files:**
- Modify: `scripts/export_pn2d_minimal6_states.py`
- Modify: `tests/regression/test_pn2d_minimal6_state_export.py`

**Interfaces:**
- Consumes: existing `prepare_exports()`, `run_exports()`, and field contract.
- Produces: manifest-bound `expected_matrix`, generalized `validate_state_matrix(states, expected_matrix)`, and captured `sentaurus_version` provenance while preserving the default six-state behavior.

- [ ] **Step 1: Add failing tests for a forty-state declared matrix**

```python
expected = tuple((topology, float(-bias))
                 for topology in ("sketch", "mirror") for bias in range(1, 21))
states = [{"topology_id": topology, "requested_bias_V": bias,
           "actual_bias_V": bias, "status": "passed"}
          for topology, bias in expected]
self.assertEqual(validate_state_matrix(states, expected), list(expected))

states[-1]["actual_bias_V"] += 2e-12
with self.assertRaisesRegex(ValueError, "does not match requested"):
    validate_state_matrix(states, expected)
```

Add a backward-compatibility test asserting that omitting `expected_matrix`
still requires exactly the legacy six states `(sketch,mirror) x (0,-12,-19)`.
Add a live-executor mock that returns `sentaurus_version="O-2018.06-SP2"` and
assert the final manifest records exactly one consistent version.

- [ ] **Step 2: Run the state-export test and verify RED**

Expected: FAIL because `validate_state_matrix()` accepts one positional argument.

- [ ] **Step 3: Implement the declared-matrix contract**

Change the signature and validation without weakening the legacy default:

```python
def validate_state_matrix(states, expected_matrix=None):
    required = set(expected_matrix) if expected_matrix is not None else {
        (topology, bias) for topology in REQUIRED_TOPOLOGIES for bias in REQUIRED_BIASES
    }
    matrix = []
    for state in states:
        topology = str(state.get("topology_id", ""))
        requested = float(state.get("requested_bias_V", math.nan))
        if state.get("status") != "passed":
            raise ValueError(f"state {topology} at {requested:g} V is not passed")
        actual = validate_final_bias(
            requested, float(state.get("actual_bias_V", math.nan))
        )
        key = (topology, requested)
        if key in matrix:
            raise ValueError(f"duplicate state {topology} at {requested:g} V")
        matrix.append(key)
        if actual != actual:
            raise ValueError("unreachable non-finite bias")
    if set(matrix) != required or len(matrix) != len(required):
        missing = sorted(required - set(matrix))
        extra = sorted(set(matrix) - required)
        raise ValueError(
            f"exact declared state matrix mismatch; missing={missing}, extra={extra}"
        )
    return matrix
```

In `prepare_exports()`, write a sorted `expected_matrix` list derived from the
requested topology and bias arguments. Thread that recorded matrix through
`run_exports()`, `validate_recovered_archive()`, and sealing/recovery call sites;
reject changes to it after preparation. Collect the executor's
`sentaurus_version` and reject mixed versions. Keep all current member hashes,
decks, TDRs, logs, imports, and failure manifests.

- [ ] **Step 4: Run state-export, Sentaurus-gate, and input tests**

Expected: all PASS, including legacy fixed-six tests and the new forty-state fixture.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/export_pn2d_minimal6_states.py tests/regression/test_pn2d_minimal6_state_export.py
D:\msys64\usr\bin\git.exe commit -m "Support declared Minimal6 Sentaurus state matrices"
```

### Task 4: Recover potential gradients and electric-field vectors

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_fields.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_fields.py`

**Interfaces:**
- Consumes: canonical SI observations and topology mesh from Task 2.
- Produces: `triangle_gradient()`, `cell_to_node_vectors()`, `cell_to_edge_vectors()`, `edge_scalar_difference()`, `mirror_vector()`, `vector_error()`, and `evaluate_field_candidates()`.

- [ ] **Step 1: Write analytic triangle, edge, mirror, zero-direction, and threshold tests**

Use the triangle `(0,0), (2,0), (0,1)` with scalar `f=3*x-4*y+2` and assert
`triangle_gradient(points, values) == (3,-4)`. Assert `E=-grad(psi)`, an x-mirror changes
`(Ex,Ey)` to `(-Ex,Ey)`, an edge projection preserves its declared orientation,
and a zero reference vector returns `DIRECTION_UNDEFINED` instead of angle zero.

- [ ] **Step 2: Run the field test and verify RED**

Expected: FAIL because `inverse_fields` does not exist.

- [ ] **Step 3: Implement exact geometry and metrics**

Use the determinant formula:

```python
def triangle_gradient(points, values):
    (x0, y0), (x1, y1), (x2, y2) = points
    f0, f1, f2 = values
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if abs(det) <= 1e-300:
        raise ValueError("degenerate triangle")
    gx = ((f1 - f0) * (y2 - y0) - (f2 - f0) * (y1 - y0)) / det
    gy = ((x1 - x0) * (f2 - f0) - (x2 - x0) * (f1 - f0)) / det
    return gx, gy
```

Area-weight cell vectors at nodes and at adjacent edges. Compare magnitude by
relative error only when the reference exceeds the declared floor; compute
direction with a clamped dot-product cosine. Emit candidates
`triangle_minus_grad_psi`, `node_area_weighted_minus_grad_psi`,
`edge_area_weighted_minus_grad_psi`, and `signed_edge_minus_delta_psi_over_h`.

- [ ] **Step 4: Run field, topology, and existing diagnostic-physics tests**

Expected: PASS with exact analytic gradients and mirror invariance.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_fields.py tests/regression/test_pn2d_minimal6_inverse_fields.py
D:\msys64\usr\bin\git.exe commit -m "Add Minimal6 field inverse candidates"
```

### Task 5: Recover quasi-Fermi-gradient and current-density semantics

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_transport.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_transport.py`

**Interfaces:**
- Consumes: Task 2 observations and Task 4 gradients/projections.
- Produces: `qf_current_density()`, `current_inverted_qf_gradient()`, `project_vector_to_edge()`, `reconstruct_edge_vector()`, `evaluate_transport_candidates()`, and typed confounding records.

- [ ] **Step 1: Write failing sign, flat-QF, inverse, projection, and confounding tests**

Assert the carrier convention implied by the repository state relations:

```python
self.assertEqual(qf_current_density("electron", 2.0, 3.0, (4.0, -5.0), q=1.0),
                 (-24.0, 30.0))
self.assertEqual(qf_current_density("hole", 2.0, 3.0, (4.0, -5.0), q=1.0),
                 (24.0, -30.0))
self.assertEqual(qf_current_density("electron", 2.0, 3.0, (0.0, 0.0), q=1.0),
                 (0.0, 0.0))
self.assertEqual(current_inverted_qf_gradient("electron", 2.0, 3.0,
                                               (-24.0, 30.0), q=1.0),
                 (4.0, -5.0))
```

When mobility is absent, assert classification `confounded` and report the
observable `mu_times_grad_qf`; when density is below floor, assert no division
occurs.

- [ ] **Step 2: Run the transport test and verify RED**

Expected: FAIL because `inverse_transport` does not exist.

- [ ] **Step 3: Implement transport candidates and support semantics**

Implement conventional-current signs explicitly:

```python
def qf_current_density(carrier, density, mobility, gradient, *, q=1.602176634e-19):
    if density < 0.0 or mobility < 0.0:
        raise ValueError("density and mobility must be non-negative")
    sign = -1.0 if carrier == "electron" else 1.0
    return tuple(sign * q * mobility * density * component for component in gradient)
```

Evaluate triangular, edge-difference, area-weighted node, area-weighted edge,
and current-inverted gradients separately for electrons and holes. Compare
Sentaurus nodal current vectors, their signed edge projections, Vela SG edge
flux, and reconstructed node/cell vectors without treating magnitudes as signed
fluxes. Report drift/diffusion terms only when their independently required
state and mobility inputs exist.

- [ ] **Step 4: Run transport, SG-flux, formula-difference, and field tests**

Expected: all PASS; flat quasi-Fermi potentials yield zero reconstructed current.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_transport.py tests/regression/test_pn2d_minimal6_inverse_transport.py
D:\msys64\usr\bin\git.exe commit -m "Add Minimal6 transport inverse candidates"
```

### Task 6: Invert avalanche driving force and reconstruct generation

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_avalanche.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_avalanche.py`

**Interfaces:**
- Consumes: Task 4 field vectors, Task 5 quasi-Fermi/current vectors, existing `van_overstraeten_alpha()`, and existing Minimal6 support integration helpers.
- Produces: `invert_van_overstraeten_alpha()`, `current_aligned_magnitude()`, `impact_generation()`, `evaluate_avalanche_candidates()`, and local/integrated generation records.

- [ ] **Step 1: Write failing forward/inverse, branch, projection, unit, and source tests**

For each electron/hole low/high branch, generate alpha with the existing forward
law and assert inversion recovers field to `1e-12` relative error. Assert alpha
outside `(0, gamma*a)` is typed below-floor or branch-ambiguous, and assert

```python
impact_generation(2.0, (3.0, 4.0), 7.0, (0.0, 24.0), q=2.0) == 89.0
```

because `(2*5 + 7*24)/2 = 89`. Add a two-triangle fixture that proves nodal
volumetric generation, triangle integration, centimetre/metre conversion, and
one-centimetre out-of-plane normalization are not mixed.

- [ ] **Step 2: Run the avalanche test and verify RED**

Expected: FAIL because `inverse_avalanche` does not exist.

- [ ] **Step 3: Implement inverse alpha and source layers**

Use the exact inverse on a known branch:

```python
def invert_van_overstraeten_alpha(alpha, *, prefactor, critical_field, gamma):
    ceiling = gamma * prefactor
    if alpha <= 0.0:
        return None, "below_numerical_floor"
    if alpha >= ceiling:
        return None, "coefficient_branch_ambiguous"
    field = -gamma * critical_field / math.log(alpha / ceiling)
    return field, "valid"
```

Reject a recovered low-branch field above the switch or a high-branch field
below it. Compare `|E|`, `|grad(phi_qf)|`, `|E dot J_hat|`, and
`|grad(phi_qf) dot J_hat|`. When electron/hole reference density is
independently declared, also evaluate `w=n/(n+n_ref)` (or `p/(p+p_ref)`) and
the interpolated driver `w*F_qf + (1-w)*|E|` without fitting the reference
density. Reconstruct native nodal generation, alpha-current generation,
triangle-integrated generation, Vela edge-partial-volume source, and node
mapping as distinct quantities with explicit units.

- [ ] **Step 4: Run avalanche, impact-ionization, cell-reconstruction, and diagnostic-physics tests**

Expected: all PASS; forward/inverse alpha agrees on both configured branches.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_avalanche.py tests/regression/test_pn2d_minimal6_inverse_avalanche.py
D:\msys64\usr\bin\git.exe commit -m "Add Minimal6 avalanche inverse candidates"
```

### Task 7: Add causal replacement and identifiability classification

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_replacements.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_replacements.py`

**Interfaces:**
- Consumes: candidate records from Tasks 4-6 and Task 1 thresholds.
- Produces: `INVERSE_DEPENDENCIES`, `run_replacement_matrix()`, `metric_summary()`,
  `classify_candidate()`, `rank_candidates()`, and `run_state_localization_control()`.

- [ ] **Step 1: Write failing dependency, closure, reverse, split, and classification tests**

Use this exact dependency order:

```python
expected = (
    "gradient_recovery", "mobility", "current_semantics",
    "impact_driving_field", "alpha_law", "geometric_integration",
    "source_to_node_mapping",
)
self.assertEqual(tuple(INVERSE_DEPENDENCIES), expected)
```

Create an analytic multiplicative fixture where every factor has a known
log10 contribution. Assert each single replacement changes only declared
downstream values, forward and reverse stages close within `1e-10 dex`, and
full replacement equals the direct target. Add classifications for one passing
candidate, two indistinguishable candidates, missing mobility confounding,
insufficient valid samples, and holdout failure.
Add a whole-state Sentaurus replay fixture and assert that it is labeled
`localization_control`, excluded from candidate ranking, and never classified
as an identified formula.

- [ ] **Step 2: Run the replacement test and verify RED**

Expected: FAIL because `inverse_replacements` does not exist.

- [ ] **Step 3: Implement deterministic replacement and classification**

Reuse the existing counterfactual engine only for arithmetic replay; keep the
new physical dependency names and availability reasons in this module. Compute
median and linear-interpolated p95 from sorted finite valid samples. Require all
configured thresholds on both combined and holdout metrics. Return
`consistent_nonunique` when passing candidates differ by <= `1e-10 dex` on all
valid Minimal6 samples, `confounded` when an independent factor is absent,
`insufficient_data` when valid support is absent, and `rejected` for gate
failures.
Implement `run_state_localization_control()` as a separate path that replaces
potential, carrier state, and quasi-Fermi state together, records downstream
localization only, and is excluded from the factor DAG, formula metrics, and
identifiability ranking.

- [ ] **Step 4: Run replacement, counterfactual, formula-difference, and contract tests**

Expected: PASS with zero unexplained residual in the analytic fixture.

- [ ] **Step 5: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_replacements.py tests/regression/test_pn2d_minimal6_inverse_replacements.py
D:\msys64\usr\bin\git.exe commit -m "Classify Minimal6 inverse replacements"
```

### Task 8: Build deterministic report, figures, CLI, and independent verifier

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_report.py`
- Create: `scripts/pn2d_minimal6_diagnostics/inverse_plots.py`
- Create: `scripts/diagnose_pn2d_minimal6_physics_inverse_audit.py`
- Create: `scripts/verify_pn2d_minimal6_physics_inverse_audit.py`
- Create: `tests/regression/test_pn2d_minimal6_inverse_report.py`

**Interfaces:**
- Consumes: all Tasks 1-7 interfaces.
- Produces: CLI arguments `--vela-root`, `--sentaurus-root`, `--supplemental-sentaurus-root`, `--out-dir`, and `--phase-base`; the complete report artifact contract; `verify_report(root)`.

- [ ] **Step 1: Write a failing two-run end-to-end report contract test**

Run the CLI twice on the same synthetic fixture and require byte-identical:

```python
ARTIFACTS = (
    "input_manifest.json", "observations_node.csv", "observations_edge.csv",
    "observations_cell.csv", "observations_contact.csv",
    "observations_integrated.csv", "candidate_metrics.csv",
    "candidate_classifications.json", "replacement_matrix.csv",
    "physics_inverse_audit.json", "physics_inverse_audit.md",
    "figure_manifest.json", "report_manifest.json", "verification.json",
    "package_manifest.json",
)
self.assertEqual({name: sha256(out_a / name) for name in ARTIFACTS},
                 {name: sha256(out_b / name) for name in ARTIFACTS})
self.assertTrue(verify_report(out_a)["passed"])
```

Require five deterministic figure pairs: `potential_field`, `qf_gradient`,
`current_density`, `alpha_generation`, and `replacement_matrix`. Mutate one CSV,
one JSON classification, one PNG pixel, and one input field in separate cases;
the verifier must reject each mutation.

- [ ] **Step 2: Run the report test and verify RED**

Expected: FAIL because the report CLI and verifier do not exist.

- [ ] **Step 3: Implement deterministic writers and figures**

Write JSON with `sort_keys=True`, `indent=2`, `allow_nan=False`, and trailing
newline. Write CSV with fixed columns, `.17g` finite formatting, and `\n` line
terminators. Fix Matplotlib backend to `Agg`, figure sizes, DPI, fonts, axis
order, legend order, metadata, and PDF creation date. After the report CLI has
written all analysis artifacts, write `report_manifest.json` last; it hashes
every raw input and report artifact except itself, `verification.json`, and
`package_manifest.json`.

The authoritative JSON must contain exact inputs, discovery/holdout keys,
thresholds, field inventory, sample-status counts, candidate metrics,
classifications, replacement closure, remote version, and production baseline.
The Markdown report must lead with the scientific conclusion and clearly
separate identified, nonunique, confounded, insufficient, and rejected results.

- [ ] **Step 4: Implement independent verification**

The verifier must not call report-building functions. It reloads raw inputs,
recomputes selected triangle gradients, current-inverted gradients, inverse
alpha, generation integrals, discovery/holdout membership, thresholds, and
replacement closure, then checks every output hash and figure pixel hash. It
writes deterministic `verification.json` only after all semantic checks pass,
then writes `package_manifest.json` hashing every input/report artifact plus
`report_manifest.json` and `verification.json`, excluding only itself. A second
verification pass validates the package manifest without changing its bytes.

- [ ] **Step 5: Run report tests and the full Minimal6 Python suite**

Run:

```powershell
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" -m unittest discover -s tests\regression -p "test_pn2d_minimal6_*.py" -v
```

Expected: all Minimal6 modules PASS and repeated report hashes match exactly.

- [ ] **Step 6: Commit**

```powershell
D:\msys64\usr\bin\git.exe add scripts/pn2d_minimal6_diagnostics/inverse_report.py scripts/pn2d_minimal6_diagnostics/inverse_plots.py scripts/diagnose_pn2d_minimal6_physics_inverse_audit.py scripts/verify_pn2d_minimal6_physics_inverse_audit.py tests/regression/test_pn2d_minimal6_inverse_report.py
D:\msys64\usr\bin\git.exe commit -m "Build Minimal6 physics inverse audit report"
```

### Task 9: Regenerate supplemental fields, execute the audit, and seal the report

**Files:**
- Modify: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`
- Generated, ignored: `build-release/pn2d-minimal6-inverse-sentaurus-fields-20260717-a/`
- Generated, ignored: `build-release/pn2d-minimal6-physics-inverse-audit-20260717-a/`
- Generated, ignored: `build-release/pn2d-minimal6-physics-inverse-audit-20260717-b/`

**Interfaces:**
- Consumes: all previous tasks and sealed Task 8 roots.
- Produces: final supplemental Sentaurus manifest, two identical report roots, independent verification, validation-doc conclusion, and a production-code no-change proof.

- [ ] **Step 1: Run all focused tests before external execution**

Run the seven new inverse-audit test modules, existing state export, diagnostic
physics, formula difference, diagnostic sweep, sweep comparison, and diagnostic
contracts modules.

Expected: every test PASS; do not start the remote run on a failing tree.

- [ ] **Step 2: Regenerate the exact forty-state supplemental Sentaurus matrix**

Run:

```powershell
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\export_pn2d_minimal6_states.py --topologies=sketch,mirror --biases=-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20 --run-id=minimal6_inverse_fields_20260717_a --output-dir=build-release\pn2d-minimal6-inverse-sentaurus-fields-20260717-a --ssh-target=sentaurus --remote-root=/home/tcad/codex_pn2d_minimal6_inverse_fields_20260717_a --importer=build-release\sentaurus_import.exe
```

Expected: forty passed states, exact requested/actual biases, one Sentaurus
version `O-2018.06-SP2`, all required QF/current/mobility/velocity/alpha/source
fields, no failed state, and complete member hashes. Preserve partial failure
manifests if the run fails; fix the cause and use a new run ID rather than
editing returned artifacts.

- [ ] **Step 3: Generate two independent report roots**

Run the report CLI twice with:

```powershell
--vela-root=build-release\pn2d-minimal6-task8-vela-fresh-20260717-a
--sentaurus-root=build-release\pn2d-minimal6-task8-sentaurus-fresh-20260717-a
--supplemental-sentaurus-root=build-release\pn2d-minimal6-inverse-sentaurus-fields-20260717-a\minimal6_inverse_fields_20260717_a
--phase-base=a5524cf
```

Use output roots ending in `-a` and `-b`. Expected: both CLIs exit zero and all
machine-readable artifacts and figure pairs are byte-identical.

- [ ] **Step 4: Run the independent verifier against both roots**

```powershell
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\verify_pn2d_minimal6_physics_inverse_audit.py --report-root=build-release\pn2d-minimal6-physics-inverse-audit-20260717-a
& "C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\verify_pn2d_minimal6_physics_inverse_audit.py --report-root=build-release\pn2d-minimal6-physics-inverse-audit-20260717-b
```

Expected: both report `passed: true`, zero hash failures, zero semantic replay
failures, mirror invariance PASS, and replacement closure <= `1e-10 dex` where
finite and nonzero.

- [ ] **Step 5: Record the evidence-backed conclusion in validation docs**

Append an inverse-audit section naming the authoritative `-a` root, input and
report hashes, exact sample counts, classifications, best formulas and their
validity domains, rejected/confounded candidates, missing evidence, and the
explicit statement that no production formula was changed. Do not describe a
candidate as identified unless it passes discovery, holdout, symmetry,
direction, magnitude, and replacement gates.

- [ ] **Step 6: Run complete regression and production-diff verification**

Run the complete Minimal6 Python suite, build the existing C++ targets, and run:

```powershell
build-release\test_fixed_state_operator_audit.exe
build-release\test_impact_ionization.exe
build-release\test_cell_reconstructed_avalanche.exe
D:\msys64\usr\bin\git.exe diff --exit-code a5524cf -- include src
D:\msys64\usr\bin\git.exe diff --check
```

Expected: all tests PASS; the `include src` diff is empty; `git diff --check`
prints nothing.

- [ ] **Step 7: Commit the validation record**

```powershell
D:\msys64\usr\bin\git.exe add docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md
D:\msys64\usr\bin\git.exe commit -m "Record Minimal6 physics inverse audit evidence"
```

- [ ] **Step 8: Request independent scientific and code review**

Provide the reviewer the design spec, this plan, commit range, authoritative
report root, supplemental root, verifier command, and explicit instruction to
classify findings as Critical, Important, or Minor. Require independent replay
of input hashes, discovery/holdout separation, formula signs and units, alpha
inversion branches, source integration, mirror transforms, figure pixels,
determinism, and the production-code no-change proof before declaring the phase
complete.
