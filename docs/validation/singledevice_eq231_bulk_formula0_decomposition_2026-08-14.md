# SingleDevice Eq. 231 bulk Formula-0 row decomposition (2026-08-14)

## Outcome

The ordinary-Silicon rows at nodes 2101, 2100, and 2044 were decomposed at
the linear and saturation endpoints.  The Formula-0 fitted nonlinear flux is
not the source of the remaining error: after a row-local unit conversion from
the independent theta=0 off-diagonal Jacobian, Vela and Sentaurus agree on
the nonlinear correction to better than `6e-7` in all three linear-endpoint
rows.

The discrepancy is the reaction closure and has two independently measured
components:

1. Vela reconstructs Lambda from the potential-like restart with a band drive
   that is about `20.7--22.2 mV` too large.  The direct Sentaurus
   `eQuantumPotential` restart value is therefore lost even at ordinary bulk
   Silicon nodes.
2. Sentaurus' reaction diagonal is exactly `1.02716991159` times Vela's at
   all three rows.  With the current configured mass ratio
   `1.06180161716`, the equivalent Eq. 231 coefficient mass is
   `1.09065067323`.

Using the direct restart Lambda and the independently recovered reaction
diagonal reduces the reconstructed full-row residual below `6e-7` at both
endpoints, without fitting any node-specific offset.

## Method

- Perturbed the union of the three rows' 13 adjacent columns by `1e-5 V`.
- Ran one full Formula-0 and one theta=0 Sentaurus NewtonPlot probe per
  column in an isolated VM directory.
- Recovered each row's unit scale from theta=0 off-diagonal derivatives only.
  Those columns contain the P1 Laplacian and do not contain the reaction
  diagonal.
- Compared the full Formula-0 derivatives and residual correction only after
  fixing that independent scale.
- Replayed every adjacent triangle separately from Vela's cell diagnostics.

The three row-local scale factors are identical to numerical precision:
`-295.290768307`, where the sign is only the opposite Sentaurus/Vela residual
convention.

## Linear endpoint

| Node | Vela full residual | Vela nonlinear correction | Sentaurus nonlinear correction | restart Lambda minus reconstructed Lambda | Sentaurus/Vela reaction diagonal | corrected full residual |
|---:|---:|---:|---:|---:|---:|---:|
| 2101 | 2.470188174 | 3.554957424 | 3.554957943 | -0.021889860 V | 1.027169911593 | -5.19e-7 |
| 2100 | 2.399042732 | 3.953150185 | 3.953150752 | -0.021185360 V | 1.027169911589 | -5.67e-7 |
| 2044 | 2.354492471 | 4.195327688 | 4.195328284 | -0.020742966 V | 1.027169911595 | -5.96e-7 |

The theta=0 non-diagonal Jacobian entries match the P1 geometry operator
after the independent scale conversion.  The Formula-0 off-diagonal entries
also match; the remaining diagonal difference is precisely the reaction
coefficient ratio above.  Per-cell contributions are individually large and
cancel across the row, so no single triangle or local edge is defective.

## Saturation endpoint cross-state check

The linear probe recovers a geometry/material coefficient, so the same single
reaction-diagonal factor was applied to the independent saturation restart.

| Node | Original Vela full residual | restart Lambda minus reconstructed Lambda | corrected full residual |
|---:|---:|---:|---:|
| 2101 | 2.496661746 | -0.022185875 V | -4.96e-7 |
| 2100 | 2.419791728 | -0.021438448 V | -5.46e-7 |
| 2044 | 2.367496318 | -0.020932337 V | -5.76e-7 |

This cross-state result rejects a one-state or one-node fit.

## Root cause in the current implementation

The potential-like restart is correct and retains the direct Sentaurus
Lambda.  During Eq. 231 assembly, however, `nodeOutputShift` and each
material-side reaction trace recompute OldSlotboom narrowing with
`deltaEg(total_impurity, n, p)`.  The generic Vela OldSlotboom implementation
uses `max(total_impurity, n, p)`.  At these inversion-layer nodes the electron
density exceeds the impurity concentration, producing about `0.091--0.093 eV`
of narrowing instead of the Sentaurus exported `0.0495621 eV`.  Half of this
difference is the observed `20--22 mV` Lambda reconstruction error.

The second issue is parameter separation.  The current
`electron_quantum_dos_mass_ratio` is used both for the DOS logarithmic drive
and for `gamma*hbar^2/(6*m*q)`.  The Jacobian probe shows that the Sentaurus
Eq. 231 reaction coefficient requires an independently configurable quantum
coefficient mass (`1.09065067323` here), while the DOS mass remains the value
used to construct the potential-like restart.

## Decision and next implementation gate

The next code change should be narrowly scoped:

1. Evaluate the Eq. 231 band-drive OldSlotboom term from total impurity only;
   do not change the already validated generic DD/BGN call sites in this
   task.
2. Separate the Eq. 231 coefficient mass from the DOS mass and retain neutral
   defaults for existing configurations.
3. Re-run this three-row oracle first.  Require a maximum absolute fixed-state
   residual below `1e-5` at both endpoints.
4. Only then rerun the two self-consistent endpoint solves; keep the 21-point
   Id-Vg curves gated until both endpoints pass.

## Artifacts

- Analyzer: `scripts/decompose_singledevice_eq231_bulk_rows.py`
- Linear JSON: `reports/eq231_bulk_formula0_decomposition_20260814/lin_bulk_formula0_decomposition.json`
- Saturation JSON: `reports/eq231_bulk_formula0_decomposition_20260814/sat_bulk_formula0_decomposition.json`
- Per-cell CSV files: `lin_bulk_formula0_cells.csv` and
  `sat_bulk_formula0_cells.csv` in the same report directory.

The licensed probe TDRs and generated exports remain under ignored build
directories and are not intended for source control.
