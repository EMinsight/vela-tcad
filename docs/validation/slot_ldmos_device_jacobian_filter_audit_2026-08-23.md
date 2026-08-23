# SLOT-LDMOS device Jacobian and residual-filter audit (2026-08-23)

## Scope

This audit followed the failed direct-bordered current-control trial at the
corrected BV branch. It intentionally did not raise the inexact-Newton forcing
cap. The checks separated the device residual/Jacobian, terminal-current row,
linear solver, continuation predictor, and residual-filter globalization.

## Root causes and fixes

| Finding | Evidence | Fix |
|---|---|---|
| Triangle-GSS residual evaluated oxide cells while local-AD and finite-difference Jacobians skipped them | Original hotspot node 6445 is shared by Si and SiO2. With impact disabled the JVP error vanished. With impact enabled, the analytic source derivative was identical to the impact-off result | `triangleGssAvalancheSourceRecordsForCell` now returns no records for nontransport materials, so primal residual and source Jacobian have identical support |
| A converged secant state could be extrapolated far off the device manifold | The old predictor changed the device residual from about `5e-6` to `23` before Newton | Predictor residual gate falls back to the last accepted state when an already-converged parent extrapolates outside the configured PDE envelope |
| Residual-filter envelopes could ratchet upward after every accepted load improvement | The failed trace accepted approximately `1.10e-6 -> 2.19e-6 -> ... -> 7.19e-5` while the load row improved only slightly | Device and load envelopes are anchored to the coupled-solve entry scales; accepted iterations cannot enlarge either reference envelope |
| Zeroing the device right-hand side used the eventual forcing cap instead of the active tolerance | A trial outside the current PDE tolerance continued to receive a tangent-only update | Zeroing is now allowed only below `coupledEffectiveEquationTolerance`; leaving that tolerance restores the normal PDE correction |
| `direct_bordered_qr` did not actually prefer QR | QR was only attempted after LU failure | The QR mode now tries rank-revealing sparse QR first, with fallback only on factorization failure |
| Current-control targets below `1e-12` were silently treated as duplicate boundary values | A requested `1e-14 A/um` increment produced only one output point | When `coupled_voltage_coefficient=0`, boundary equality uses the configured locator tolerance rather than the voltage-sweep hard-coded tolerance |

## Jacobian evidence

1. At old hotspot node 6445, after excluding nontransport cells from the
   Triangle-GSS residual, the physical-potential JVP error at a `1e-6 V`
   perturbation dropped from about `4.414e-1` to `1.678e-10`.
2. At corrected hotspot node 10462, all 15 basis directions in the one-cell
   vertex star (five nodes times three potential blocks) matched central finite
   differences. The maximum relative error at `1e-6 V` was `3.282e-7`.
3. The direct bordered circuit row is internally consistent: representative
   analytic versus finite-difference terminal-current directional derivatives
   differed by only a few parts in `1e5` to `1e4`, while the augmented linear
   residual remained near `1e-27` to `1e-11` depending on row scaling.
4. LU and forced rank-revealing QR produced the same physical tangent, ruling
   out the linear backend as the source of the failed nonlinear trial.

The full near-null tangent JVP is not a stable relative-error metric here: its
analytic device action is deliberately close to zero, so residual subtraction
noise dominates the finite-difference quotient. The resolved local basis
audits above are the deterministic matrix-entry qualification.

## Corrected branch evidence

The stale high-voltage state was reclosed at fixed drain voltage after the
material-support fix:

| Drain voltage (V) | Drain current (A/um) | Status |
|---:|---:|---|
| 15.856737161516595 | 4.244286026391734e-10 | converged |
| 15.857737161516594 | 4.244433552675878e-10 | converged |

The measured local slope is approximately
`1.47526284144e-11 A/(um V)`. The earlier locator target around
`3.3236e-9 A/um` came from a checkpoint containing the oxide-side pseudo
avalanche source and must not be reused.

## Current status

The corrected fixed-voltage branch and local Jacobian are qualified. The
residual filter now rejects or bounds bad collective directions instead of
ratcheting the PDE residual. A bounded direct-current smoke probe still does
not close its current row within eight iterations: on this very flat branch,
the avalanche residual has strong collective-direction curvature and the
accepted step is limited by the solve-entry PDE envelope.

Therefore the previous IALMob A/B high-voltage checkpoints and BVDS result are
invalidated by the primal source-support correction. Before extracting BVDS,
both A/B branches must be regenerated from a low-voltage state with the
corrected source. If direct current control remains necessary after that
rebuild, its next solver task is an avalanche-driving-field-aware trust region
or an explicit device-manifold corrector, not a larger forcing cap.

## Verification

| Check | Result |
|---|---|
| `test_coupled_load_line` | 65 assertions, 9 cases passed |
| `test_dc_sweep` | 3471 assertions, 98 cases passed |
| `test_mos_mixed_material` | 1439 assertions, 10 cases passed |
| `test_newton_solver "[triangle_gss]"` | 8 assertions, 2 cases passed |
| `slot_ldmos_jvp_audit_preparation` | passed |
| Corrected fixed-voltage branch probe | 2/2 points converged |
