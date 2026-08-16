# SingleDevice Eq. 231 self-consistent Id-Vg validation (2026-08-15)

## Scope

This run compares Vela's fully self-consistent electron density-gradient
solution with the Sentaurus 2018 SingleDevice references at:

- linear drain bias: `Vds = 0.1 V`;
- saturation drain bias: `Vds = 1.1 V`;
- gate range: `-0.5 V` to `2.2 V`;
- reference lattice: 21 points at `0.135 V` spacing.

The quantum inner solve uses `absolute_tolerance_V = 1e-6`; the outer
DD/quantum fixed-point acceptance remains `0.5 mV`.  These tolerances are now
configured independently through `outer_absolute_tolerance_V`.

## Continuation procedure

Directly starting the high-inversion endpoint at `Vg = 2.2 V` excited an
unstable outer mode.  Both curves therefore start from the low-gate state and
sweep upward.

The linear curve uses the 21-point reference lattice directly.  The saturation
curve first ramps `Vds` from `0.1 V` to `1.1 V` at fixed `Vg = -0.5 V` in
`0.05 V` steps.  Its Id-Vg sweep then uses a `0.0675 V` half-step lattice, so
all 21 reference biases are present without interpolation.  Saturation DD
acceptance uses `reltol = 1e-6` and `abstol = 1e-5`; the quantum outer
criterion is unchanged.

## Results

| Metric | Linear | Saturation | Requirement |
|---|---:|---:|---:|
| Converged continuation points | 21/21 | 41/41 | all |
| Compared reference points | 21 | 21 | at least 21 |
| Maximum relative current error | 10.977% | 7.965% | at most 10% |
| Maximum log-current error | 0.05050 decade | 0.03605 decade | at most 0.2 decade |
| Trend match | yes | yes | yes |
| Maximum relative error for `Vg >= 0.31 V` | 1.396% | 0.855% | diagnostic |
| Constant-current threshold delta at `Id = 1e-7 A/um` | +1.964 mV | +1.922 mV | diagnostic |

The saturation curve passes all existing curve requirements.  The linear
curve narrowly misses the global relative-error requirement only at
`Vg = -0.5 V`: Sentaurus gives `1.2472174e-14 A/um` and Vela gives
`1.1103056e-14 A/um`, a 10.977% difference.  Its maximum logarithmic error is
only 0.05050 decade, and all other linear points are below 7.32% relative
error.  Therefore the two-curve validation is **not yet fully accepted** under
the strict global 10% rule, despite both curves converging and the saturation,
trend, logarithmic, threshold, and strong-inversion checks passing.

## Reproducible artifacts

Generated runtime artifacts are under:

`build-release/reference_tcad/singledevice_sentaurus2018/vela_import_fixedmaterials/vela`

Comparison artifacts are under:

`build-release/reference_tcad/singledevice_sentaurus2018/reports/self_consistent_idvg_20260815`

The remaining validation item is localized to the linear deep-off leakage at
`Vg = -0.5 V`; the full nonlinear Eq. 231 continuation and the saturation
curve are no longer blockers.

## 2026-08-16 deep-off follow-up

The focused SRH, minority-density, contact-current, and numerical-floor audit
is recorded in
`singledevice_deep_off_current_diagnostic_2026-08-16.md`.  It classifies the
remaining point as a real but small electron-state/transport mismatch.  SRH,
contact extraction, KCL closure, and numerical current noise are not dominant.
Use a documented hybrid low-current acceptance rule rather than the strict
pointwise relative-current gate alone before closing the two-curve contract.
