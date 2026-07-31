# PN2D BV M0/M2 template-default acceptance v2

Date: 2026-07-31

Status: complete; IIC convergence fixed; SG/Laux default proposal rejected and
rolled back.

## Decision

The v2 machine decision is:

```text
pn2d_bv_template_default_not_accepted
```

`production_default_change_authorized=false`.  SG/Laux remains available as
the named opt-in profile, while the PN2D BV template default is restored to
`legacy_cell_reconstructed`.  The generated evidence is under:

```text
build-release/pn2d-bv-template-default-prospective-v2-default-20260731/
```

The unified decision is `acceptance/acceptance.json`.

## IIC carrier-row convergence correction

The previous M0/M2 IIC failures occurred after the actual Newton residual had
already reached about `1e-12`.  The carrier-row acceptance diagnostic used the
observable impact-ionization source even when `coupling_mode=postprocess_only`.
That source was not present in the solved continuity residual or Jacobian, so
the diagnostic rejected an otherwise valid avalanche-off-equivalent state.

The assembler now exposes separate semantics:

- `carrierContinuityTermDiagnostics`: observable terms, including the IIC
  source;
- `carrierContinuityEquationTermDiagnostics`: only terms present in the
  equations solved by Newton.

Newton row scaling, carrier-row convergence, global closure acceptance, and
their trace output use the equation-term interface.  Postprocess-only source
export remains unchanged, and self-consistent avalanche continues to include
the source in both the equation and diagnostics.

Direct repair validation completed all 29 exact points through `-20 V` for
both M0 and M2.  The repaired IIC IV and state hashes equal avalanche-off at
every requested bias, as required by the postprocess-only semantics.

## Frozen v2 policy

The v2 contract was frozen before the accepted default-render simulations:

```text
docs/validation/contracts/pn2d_bv_m0_m2_template_default_acceptance_v2.json
SHA-256 79ccf6f6427b1fda4890581fe10299fcd9d6f8f3f8d986527936129702c275d2
```

All numerical curve, gain, knee, closure, exact-lattice, and determinism
thresholds are unchanged.  The only estimator-domain clarification is
prospective:

- both simulators have a slope crossing: compare `V_slope` with the existing
  `0.10 V` threshold;
- neither has a crossing in the frozen window: accept the typed outcome
  `shared_no_slope_crossing_in_frozen_window`;
- only one has a crossing: fail closed.

The evaluator also rejects a profile override as default-value evidence.  The
final runs therefore used a temporary atomic SG/Laux template default and a
fresh render whose manifest contains no profile override.  The template was
rolled back after the failed decision; its restored SHA-256 is
`c525e0b78d3142bdf37086dbb70086ead54a48f3a276000c25286294bc9b372b`.

## Results

Both M0 and M2, in both independent runs, completed avalanche-off,
IIC/postprocess-only, and avalanche-on on all 29 requested points.  IV,
physics configuration, process-probe, and all state hashes are deterministic.

| Metric | M0 | M2 | Limit | Result |
|---|---:|---:|---:|---|
| effective curve median error (dex) | 0.002030 | 0.046355 | 0.05 | pass / pass |
| effective curve P95 error (dex) | 0.003708 | 0.075860 | 0.10 | pass / pass |
| effective curve maximum error (dex) | 0.003749 | 0.078827 | 0.15 | pass / pass |
| effective gain median error (dex) | 0.002186 | 0.041009 | 0.05 | pass / pass |
| effective gain maximum error (dex) | 0.003702 | 0.073480 | 0.10 | pass / pass |
| knee median error (dex) | 0.003045 | 0.057983 | 0.05 | pass / **fail** |
| knee maximum error (dex) | 0.003749 | 0.078827 | 0.10 | pass / pass |
| absolute `V_break` error (V) | 0.001 | 0.011 | 0.10 | pass / pass |
| `V_slope` outcome | both present, 0.000819 V error | shared no crossing | policy | pass / pass |
| adjacent-slope RMSE (dex/V) | 0.002659 | 0.055059 | 0.20 | pass / pass |

The maximum avalanche-on global continuity closure ratios also pass:

| Level | Electron | Hole | Limit |
|---|---:|---:|---:|
| M0 | 7.13696e-6 | 7.12422e-6 | 0.01 |
| M2 | 4.70138e-7 | 2.15322e-6 | 0.01 |

## Interpretation and next work

The IIC carrier-row convergence issue is closed and is no longer a blocker.
The shared absence of an M2 slope crossing is also no longer an ambiguous
failure.  The remaining default-value blocker is narrower: the M2 same-grid
Sentaurus comparison has a knee-domain median log-current error of
`0.0579825 dex`, exceeding the frozen `0.05 dex` limit by `0.0079825 dex`.

The next investigation should therefore be observation-only and M2-local:
rank the bias-resolved error onset from `-18 V` to `-20 V` against the frozen
SG/Laux source components, edge currents, driving fields, mobility, and source
mapping.  No default change or threshold relaxation is authorized until a new
prospective contract is frozen and all gates pass.
