# PN2D BV Current Progress Summary

Prepared for external model-assisted analysis.

## Current Goal

The PN2D Sentaurus2018 BV comparison has moved past the basic continuation
question. Vela can now complete the `-20 V` sweep under the stabilized
continuation setup, but the curve still does not fully reproduce the Sentaurus
terminal avalanche knee.

Sentaurus reference behavior in the `-20..-10 V` window:

- 1 V current growth ratio `>1.5` appears at `-19 V`.
- 1 V current growth ratio `>2.0` appears at `-20 V`.

## Key Findings

### Low-Bias Continuation Blocker

The low-bias inability to continue was not fixed by Charon-style minimum-field
or current-aligned driving-force knobs:

- `minimum_field_V_m=5.0e6` did not recover the `-3 V` gate.
- `driving_force=effective_field_parallel_j` did not recover the `-3 V` gate.
- Branch guard detects a real electron-density branch jump, but is diagnostic
  rather than a continuation strategy.

The immediate low-bias blocker was a quasi-Fermi update/branch interaction.

### QF Update Cap

`quasi_fermi_update_limit_V=0.05` is the current stabilizing continuation
parameter.

| `quasi_fermi_update_limit_V` | outcome |
|---:|---|
| `0.1` | fails near `-0.3268359 V` with `electron_density_p95_jump_exceeded` |
| `0.05` | reaches `-20 V`, `491/491` converged rows |
| `0.025` | too restrictive; low-bias `line_search_non_decrease` |

However, `qf_limit=0.05` is a continuation-stability fix, not an avalanche-knee
fix. With the default source support it has max log10 current error
`0.891695` decades and no `>1.5` or `>2.0` knee marker.

## Source Ownership Scan

All rows below use:

- `quasi_fermi_update_limit_V=0.05`
- secant predictor
- carrier-density branch guard

| source setting | result | max log10 current error, `-20..-10 V` | knee marker |
|---|---:|---:|---|
| default / half support | reaches `-20 V` | `0.891695` decades | none |
| `source_volume_factor=0.875` | reaches `-20 V` | `0.453746` decades | none |
| `source_volume_factor=0.90625` | reaches `-20 V` | `0.381328` decades | none |
| `source_volume_factor=0.921875` | reaches `-20 V` | `0.349196` decades | `>1.5` at `-20 V`, no `>2.0` |
| `source_volume_factor=0.9375` | fails at low bias | n/a | `-0.6897608 V` then `line_search_non_decrease` |
| `source_volume_policy=edge_box` | fails at low bias | n/a | `-0.4148548 V` then `electron_density_p95_jump_exceeded` |
| `edge_box` with p95 guard relaxed to `2.1 dex` | still fails | n/a | `-0.4469347 V` then `line_search_non_decrease` |

Current stable scalar-source boundary is between `0.921875` and `0.9375`.
`0.921875` is the best stable scalar candidate so far, but it still does not
recover the Sentaurus `>2.0` breakdown marker.

## Excluded Likely Cause

The SG avalanche source-derivative Jacobian has already been probed:

- synthetic block relative difference is about `2e-6`
- real `-20 V` state finite-difference block checks are aligned

Current evidence does not point to the SG avalanche source Jacobian as the BV
knee blocker.

## Open Analysis Questions

1. Why does `source_volume_factor=0.9375` fail near `-0.6897626 V` with a small
   residual but a very large raw Newton step?
   - failure residual norm: `2.25599e-09`
   - step norm: `47.3865`
   - carrier densities remain positive and finite

2. Is there a stable scalar factor in `0.921875..0.9375` that recovers the
   Sentaurus `>2.0` marker?

3. Is full `edge_box` too aggressive, requiring an intermediate ownership
   geometry rather than a scalar source-volume factor alone?

4. Is a stronger continuation strategy needed?
   - current-aware predictor
   - pseudo-arclength continuation
   - low-bias branch-state regularization
   - more localized minority quasi-Fermi update caps

## Important Artifacts

- Validation document:
  `D:\code-repo\vela-tcad\docs\validation\pn2d_bv_validation.md`
- Latest high source scan:
  `D:\code-repo\vela-tcad\build-release\reference_tcad\pn2d_sentaurus2018\reports\qflim0p05_high_source_scan`
- `0.9375` failure diagnostics:
  `D:\code-repo\vela-tcad\build-release\reference_tcad\pn2d_sentaurus2018\reports\qflim0p05_high_source_scan\factor_0p9375_newton_failure_diagnostics.json`
- Initial source factor scan:
  `D:\code-repo\vela-tcad\build-release\reference_tcad\pn2d_sentaurus2018\reports\qflim0p05_source_factor_scan`

## Latest Verification

The latest checks passed:

- `git diff --check`
- `ctest --test-dir build-release --output-on-failure -R impact`, 5/5 passed

## Recommended Next Step

Binary-search the scalar factor interval `0.921875..0.9375`, starting with
`0.9296875`, while preserving the same branch guard. In parallel, inspect the
first failing Newton transition on the unstable side to determine whether the
failure is a continuation-method limitation or a source-ownership modeling
limit.

## 2026-07-14 Minimal6 Real Six-State Audit Status

Task 6 is complete on the recovered real state root
`build-release/reference_tcad/pn2d_sentaurus2018_minimal6/state_exports/minimal6_states_live_20260713_v2`.
All six deterministic production C++ replays exited 0, replay provenance was
freshly verified, and the joined report passed with `36/54/24` rows and
`7 PNG + 7 PDF` figures.

The strict results are: imported-state parity error `0 < 1e-12` and maximum
C++/Python formula error `1.7995489542954602e-12 < 5e-12`. The real report is
`build-release/reference_tcad/pn2d_sentaurus2018_minimal6/reports/minimal6_fixed_state_audit_20260714`.

This closes the minimal6 fixed-state operator-audit workflow only. No Vela
Newton, Gummel, or continuation solve was run, and the three fixed biases are
not a BV sweep. Sentaurus/Vela current and avalanche-source differences remain
order-one diagnostics with no acceptance threshold; at `-12 V` and `-19 V`,
the summed Vela sources are approximately `4.4944e14` and `3.0765e18 m^-1 s^-1`,
versus Sentaurus `3.1472e2` and `2.8566e5 m^-1 s^-1`. The next physical task
remains current/source-semantics reconciliation before any minimal6 breakdown
claim or nonlinear Vela BV validation.

## 2026-07-15 Minimal6 Tasks 1-8 completion status

Tasks 1-8 are complete for the scoped diagnostic workflow. This supersedes the
previous statement that no nonlinear Vela solve had been run.

- Six-state state/replay evidence: PASS for sketch/mirror x `0/-12/-19 V`;
  manifest SHA-256
  `43305a3b4b3f565d600c2ccb783af36c62ed385e0ea6aa9c13f71b3a34e370c0`.
- Recovery evidence: 355 members, zero hash mismatches; seal SHA-256
  `c343539775337038848adc3a1f88b45880110d6960a3fe8ea3f3309032645f6d`.
- Fixed-state formula report: exact `36/54/24` identities and five reviewed
  figure pairs; strict result remains `insufficient_data` because named
  counterfactual substitutions lack raw Vela state input.
- Vela sweep: sketch and mirror accept exact `-1 V`; both first fail at
  `-1 -> -2 V` with exit code 1 and `nonfinite_residual`. Native/reconstructed
  production SG source totals close exactly at the accepted endpoints.
- Sentaurus sweep: sketch and mirror each accept all 21 exact targets from
  `0` through `-20 V`; no failed transition. The regenerated clean 0 V imports
  pass the corrected scalar/vector field selector.
- Final comparison: the only common exact bias is `-1 V`. The terminal-current
  ratios are `1.4457180849345512e-4` (sketch) and
  `1.5145723823453054e-4` (mirror), both below the recorded leakage-like
  threshold. Forty Sentaurus-only checkpoints are retained without
  interpolation or tail extrapolation.

Primary generated roots are
`build-release/pn2d-minimal6-formula-diff-task8-20260715/`,
`build-release/pn2d-minimal6-vela-task8-final-r3-20260715/`,
`build-release/pn2d-minimal6-sentaurus-task8-final-r2-20260715/`, and
`build-release/pn2d-minimal6-comparison-task8-final-r2-20260715/`.

The physical conclusion has not changed: the evidence exposes a large
cross-solver operator/state gap but does not identify a closed causal ranking
or a physical breakdown voltage. Every curve is labelled:

> minimal6 diagnostic sweep; not a physical BV curve
