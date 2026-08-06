# BVmethods NMOS SG edge and path-ionization follow-up (2026-08-04)

## Scope

This follow-up fixes the 6.4 V state and audits, in order, the electron/hole
SG fluxes, ionization coefficients, geometric weights, source mapping, and the
Sentaurus `ComputeIonizationIntegrals` / `BreakAtIonIntegral(3 1.)` semantics.

## Edge-source ledger at 6.4 V

The matched 3970-edge ledger is written under
`btbt_e2_iic_edge_ledger_20260804`.

- Rebuilding the Vela source as `alpha * flux * edge_area` reproduces the
  stored source to roundoff (ratio `0.9999999999999997`). There is no missing
  factor of two in the edge-to-node scatter or the 2-D geometry.
- Replacing only Vela alpha by the imported Sentaurus alpha increases the
  source by `1.3733234x`.
- Replacing only the edge-projected current by the imported Sentaurus current
  increases it by `1.0253550x`.
- Using Sentaurus alpha and the Sentaurus current-vector magnitude increases
  it by `2.1104607x`; direct imported Sentaurus impact generation gives
  `2.0958311x`. Their approximately 0.7% agreement identifies current-vector
  support plus alpha distribution, rather than source mapping, as the main
  integrated-source gap.

## Vela current-magnitude control matrix at 6.4 V

All rows reuse the same accepted 6.4 V state.

| alpha driving force | current magnitude | source ratio to existing QF/edge baseline |
|---|---:|---:|
| quasi-Fermi gradient | edge scalar | 1.0000 |
| Eparallel | edge scalar | 0.2348 |
| quasi-Fermi gradient | existing dual-face vector | 3.8555 |
| Eparallel | existing dual-face vector | 3.7937 |

The unguarded dual-face reconstruction overshoots the approximately 2.10x
Sentaurus vector target and is therefore not promoted to a global default.

## Path-integral implementation

`PathIonizationIntegral` implements the Sentaurus Device User Guide O-2018.06
Eqs. 469-470 on ordered, piecewise-constant path segments. It uses stable
log-sum-exp accumulation and traces monotone-potential paths through local
electric-field maxima. Paths are ranked by the arithmetic mean of electron-
and hole-injection integrals.

The DC-sweep diagnostic is configured as:

```json
"path_ionization_integrals": {
  "enabled": true,
  "csv_file": "path_ionization_integrals.csv",
  "max_paths": 3,
  "break_rank": 3,
  "break_value": 1.0,
  "driving_force": "solver"
}
```

`break_rank: 3` stops fixed-point, arclength, and adaptive sweeps when the
third-largest mean path integral reaches one. The path driving force is
independent of the avalanche-generation source driving force, so an
approximate-breakdown diagnostic does not silently change the continuity
equations.

## Unit root cause

Vela `unit_scaling` deliberately combines micrometer coordinates with
inverse-centimeter ionization coefficients. A dimensionless alpha-length
product is therefore not the raw internal product:

`alpha[m^-1] * length[m] = alpha_internal * length_internal * 1e-4`.

The old VTK local alpha-length proxy and the first path implementation omitted
this conversion and were exactly `1e4` too large. The implementation now
converts both quantities to compatible SI units before accumulation. The
legacy local fields remain available, while new VTK fields use the explicit
`ElectronPathIonIntegral`, `HolePathIonIntegral`, and
`MeanPathIonIntegral` names.

## Driving-force comparison at 6.4 V

For the third-ranked path:

| path driving force | electron integral | hole integral | mean integral |
|---|---:|---:|---:|
| ElectricField | 1.040661 | 1.069337 | 1.054999 |
| GradQF | 0.908447 | 0.871842 | 0.890144 |
| Vela Eparallel projection | approximately 0 | 0.542118 | 0.271059 |

`ElectricField` is the only tested branch that produces an integral near one
at the independently extracted 6.377494 V current-IIC crossing. It is the
driving force the Sentaurus manual recommends for field-only approximate
breakdown. However, the official coupled `pp4` deck explicitly uses
`Avalanche(Eparallel)` with `AvalPostProcessing`; therefore the Eparallel row,
not the ElectricField row, is the like-for-like path reference for that deck.
At the highest currently accepted Vela state, 7.0 V, the Eparallel third-ranked
mean integral is `0.290136`, so the official rank-3 threshold has not been
reached anywhere on the validated 0-7 V branch.

## ElectricField rank-3 crossing (diagnostic control only)

Exact accepted Vela states give:

| bias (V) | rank-3 mean integral | seed edge |
|---:|---:|---:|
| 6.0 | 0.9765795 | 3520 |
| 6.1 | 0.9961454 | 3520 |
| 6.2 | 1.0157368 | 3520 |
| 6.3 | 1.0353548 | 3520 |
| 6.377494 | 1.0505759 | 3520 |
| 6.4 | 1.0549991 | 3520 |

Linear interpolation between 6.1 V and 6.2 V gives an ElectricField-control
Vela rank-3 crossing of `6.119675 V`.

Two Sentaurus references must remain separate:

- `6.377494 V` is the postprocessed current-IIC condition
  `q*Integral(AvalancheGeneration)=abs(Idrain)`; it is not a path threshold.
- The original `pp4_des.cmd` contains `BreakAtIonIntegral(3 1.)`, and
  `n4_des.log` terminates at the final drain point `10.448267 V` because the
  ionization integrals meet that specification.

Because this control uses ElectricField while `pp4` uses Eparallel, the
`6.119675 V` and `10.448267 V` values are not a like-for-like error metric.
They demonstrate why the path driving force must remain explicit and why
current-IIC and path-integral criteria cannot be conflated. The equations,
units, ranking, and stop semantics are implemented; numerical closure of the
coupled Eparallel path still requires extending the accepted Vela branch to
the Sentaurus path-stop range. The 6.377494 V current-IIC comparison remains a
separate source/current integration task.

## Remaining closure work

1. Extend the accepted Vela Eparallel branch beyond 7 V toward the official
   `10.448267 V` path stop, without enabling self-consistent avalanche source
   feedback.
2. Export the Sentaurus top-three path coordinates and per-path
   `eIonIntegral`, `hIonIntegral`, `eAlphaAvalanche`, and `hAlphaAvalanche`,
   then compare them with Vela seed edge 3520 path coordinates.
3. Replace graph-edge steepest descent by cell-gradient streamline tracing
   with geometric direction continuity and an explicit depletion-zone stop
   criterion if the path geometry does not align.
4. Reconcile the Vela Eparallel projection, which suppresses the electron
   coefficient on the current top paths.
5. Preserve separate acceptance gates: path rank-3 stop versus the
   `q*Integral(AvalancheGeneration)=abs(Idrain)` current-IIC crossing. Only
   after both close should the source operator be reused for resistor and
   voltage-to-current methods.

## Verification

- `test_path_ionization_integral`: 17 assertions in 4 test cases pass.
- `test_impact_ionization [impact][diagnostic][vtk]`: 112 assertions in 3
  test cases pass.
- DC step-control tests pass, including explicit diagnostic stop handling.
- The 6.4 V source-factorization, four-mode current control, and three-mode
  path-driving-force runs all converge from the same accepted state.
- The exact 7.0 V Eparallel path diagnostic converges and remains below the
  rank-3 threshold (`0.290136 < 1`).

## 10.448 V current-line repair update

The continuous-cell tracer now distinguishes three path tangents:
electrostatic field, reconstructed SG carrier current, and the recovered
quasi-Fermi gradient used only as an explicit diagnostic.  Previously the
`electron_current` name incorrectly selected the quasi-Fermi-gradient vector.
That vector is well behaved on the main conducting branch but is dominated by
roundoff at the very-low-current local field maxima, causing several false
paths to merge into the high-field drain corridor.

For SG-current tracing, nodal current vectors below `1e-8` of the global nodal
current maximum now fall back to the electrostatic-field direction. This
changes path geometry only; the seed/stop field, Eparallel driving values,
ionization coefficients, and Sentaurus Eq. 469-470 integration are unchanged.

At the official final bias `10.4482667308 V`, the repaired top-three result is:

| Rank | seed node | Vela electron | Vela hole | Vela mean | Sentaurus mean | mean error |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 990 | 1.544439 | 2.178582 | 1.861511 | 1.794564 | +3.7% |
| 2 | 1460 | 1.443188 | 1.911437 | 1.677313 | 1.544187 | +8.6% |
| 3 | 1381 | 1.366210 | 1.599848 | 1.483029 | 1.428772 | +3.8% |

The seed identities and ordering now match the three positive Sentaurus TDR
plateaus exactly. The next-highest Vela path is only `2.73e-4`, so low-current
roundoff paths no longer contaminate `BreakAtIonIntegral(3 1.)`. The remaining
4--9% mean-integral differences are coefficient/path quadrature errors, not a
missing node-1460 channel or rank-selection error.

The `1e-8` reliability floor is not a point fit: `1e-9` gives the same three
means to better than `2e-8` and keeps the fourth path at `2.71e-4`. At `1e-7`
and `1e-6` the same three seed identities remain isolated, although the node
990 path begins to change as physically resolved current vectors are replaced
by electric-field tangents. This brackets `1e-8` on the stable SG-current side
of the transition.

The ordered quadrature is now oriented along the electric field, from higher
to lower electrostatic potential. The former lower-to-higher ordering reversed
the `x` convention in Eqs. 469-470 and increased all three mean errors.

Reproducible repaired output:
`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_trace_oriented_20260804/postprocess_only`.

## Full corrected high-voltage branch recomputation

The repaired tracer was rerun uniformly on every accepted state from 7.0 V
through the official final voltage.  This run includes SG-current path
geometry, the `1e-8` current-direction reliability floor, high-to-low
potential integration, and nodal-seed corridor deduplication.  The latter
keeps one physical streamline for each nodal local maximum instead of ranking
multiple incident-cell launches from the same seed as separate paths.

- All 36 requested states converged: 7.0--10.4 V in 0.1 V increments, plus
  10.4482667308 V.
- The top-three seed set is `{990, 1460, 1381}` at all 36 states.  Node 1381
  is the third-ranked path throughout the branch.
- The node-1381 integral is monotone on the complete branch.  The largest
  fourth-ranked result is only `2.80101e-4`, so a duplicate or low-current
  side path cannot satisfy `BreakAtIonIntegral(3 1.)`.
- The first two paths exchange order near 9.0 V and again at the last two
  states, but this does not alter the rank-3 threshold.

Selected rank-3 values are:

| drain bias (V) | seed node | Vela mean integral |
|---:|---:|---:|
| 7.000 | 1381 | 0.95551895 |
| 7.200 | 1381 | 0.98754601 |
| 7.277 | 1381 | 0.99996047 |
| 7.278 | 1381 | 1.00012183 |
| 7.300 | 1381 | 1.00367218 |
| 8.000 | 1381 | 1.10476361 |
| 9.000 | 1381 | 1.25638770 |
| 10.000 | 1381 | 1.41197065 |
| 10.448267 | 1381 | 1.48302913 |

Linear interpolation of the separately resolved 7.277 V and 7.278 V states
places the current Vela arithmetic-mean rank-3 threshold at
`7.277244958 V`.  This is a smooth crossing on one persistent node-1381
path, not an interpolation across a path identity change.

The Sentaurus reference has materially different path-appearance semantics.
In the original adaptive run, the last sub-threshold-topology state is
8.595 V and the next accepted state is 10.448267 V; the latter terminates the
sweep.  A separate `WriteAll` run at approximately 0.1 V spacing shows only
two distinct high-integral corridors at 9.151 V (path numbers 0 and 1 are
duplicates of one corridor), while a third high-integral corridor appears by
9.251 V.  Thus the dense Sentaurus data do not bracket a smooth third-path
integral crossing: they bracket a topology/seed-appearance event.  The
official 10.448267 V stop is additionally an overshoot caused by its final
8.595-to-10.448 V adaptive step.

Consequently, `7.277244958 V` is the completed Vela branch result, but it is
not yet a closed Sentaurus result.  Vela preserves three high-field paths from
7.0 V onward, whereas Sentaurus exposes the third distinct corridor only near
9.2 V.  Endpoint magnitudes at 10.448267 V agree within 3.7--8.6%; the
remaining breakdown-voltage discrepancy is dominated by the voltage at which
the third path is created/accepted, not by the final integral magnitude.

Sentaurus documentation names the ranking quantity `I_mean` but does not give
an explicit algebraic definition in the approximate-breakdown section.  Vela
currently ranks the arithmetic mean of Eqs. 469 and 470.  The exported
Sentaurus `MeanIonIntegral` field is not exactly the pointwise arithmetic mean
of the exported electron and hole fields, so no unverified conversion formula
has been introduced.  The threshold above is therefore explicitly reported
as the current Vela arithmetic-mean criterion.

Reproducible outputs:

- complete branch:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_oriented_final_full_branch_20260805`;
- dense crossing:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_oriented_dense_crossing_7p2_7p3_20260805/postprocess_only`;
- dense Sentaurus `WriteAll` log:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_eparallel_vector_exact_20260804/extracted/iic_multibias_des.log`.

After the nodal-seed deduplication update, `test_path_ionization_integral`
passes 39 assertions in 8 cases, the path-ionization DC-sweep tests pass 23
assertions in 3 cases, and the impact-ionization diagnostic tests pass 327
assertions in 16 cases.
