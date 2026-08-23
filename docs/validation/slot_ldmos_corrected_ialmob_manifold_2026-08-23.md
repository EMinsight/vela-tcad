# SLOT-LDMOS corrected IALMob device-manifold continuation

Date: 2026-08-23

## Outcome

The repaired Newton globalization and previous-state device-manifold
continuation produced strict IALMob off/on pairs through 1 V.  Both cases use
the same mesh, avalanche model (`local_ad` source Jacobian), nonlinear
tolerances, block residual filter, and uniform quasi-Fermi trust region.  The
only controlled physics difference is `masetti_field` versus
`masetti_field_lombardi` at the `Silicon_1/Oxide_1` interface.

| Segment | Bias points | off | on | Rejected production points |
|---|---|---:|---:|---:|
| Dense low voltage | 0.06--0.10 V, 10 mV step | 5/5 | 5/5 | 0 |
| Post-dense | 0.12, 0.15 V | 2/2 | 2/2 | 0 |
| 0.20 V recovery | 0.16--0.20 V, 10 mV step | 5/5 | 5/5 | 0 |
| 1 V extension | 0.25, 0.30, 0.40, 0.50, 0.75, 1.00 V | 6/6 | 6/6 | 0 |

The 1 V extension internally used accepted 50 mV substeps when crossing the
larger requested landmark gaps.  Every recorded 0.16--0.20 V point converged
in four Newton iterations and had zero quasi-Fermi-bound violations.

## Paired current results

| Bias (V) | IALMob off Id (A/um) | IALMob on Id (A/um) | on/off |
|---:|---:|---:|---:|
| 0.20 | 1.791391361e-11 | 1.707829487e-11 | 0.95335 |
| 0.50 | 3.730091517e-11 | 3.450552136e-11 | 0.92506 |
| 1.00 | 6.427589711e-11 | 5.870505366e-11 | 0.91333 |

The BVDS criterion of `1e-7 A/um` has not been reached, so these results do
not yet define BVDS or the IALMob BVDS shift.

## High-voltage predictor ablation

A 1.00-to-1.25 V first jump with a bounded secant deck failed immediately
with `line_search_non_decrease` and zero accepted Newton updates.  A revised
deck starting at 1.05 V converged at 1.05 and 1.10 V, but the secant-predicted
1.20 V target generated repeated rejected attempts (4, 6, and 8).  The run was
stopped before exhausting 26 retries.

This establishes that secant state extrapolation remains unsafe for the
production avalanche branch even with an extrapolation ratio capped at one.
The prepared production high-voltage deck therefore uses:

- no state extrapolation;
- a 50 mV initial continuation step;
- growth factor 1.25;
- a 250 mV maximum step;
- unchanged nonlinear forcing and convergence tolerances.

The next run must continue both cases from the paired 1 V states to at least
12 V, then use a denser near-fold segment and the coupled series-resistor
boundary to cross the fixed-voltage fold and bracket `1e-7 A/um`.

## Verification

- Python generator and analyzer regression: 13 tests passed.
- `test_newton_solver`: 1281 assertions in 94 cases passed.
- `test_line_search_backtrack_failure`: 42 assertions in 6 cases passed.
- `test_dc_sweep`: 3471 assertions in 98 cases passed.
- `git diff --check`: passed (line-ending warnings only).
