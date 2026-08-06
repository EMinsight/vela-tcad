# BVmethods NMOS nodal-colocated path audit (2026-08-05)

## Scope

This audit reruns the Sentaurus-style path ionization diagnostic using the
same nodal physics already validated by `nodal_eparallel_p1`: complete
transport-vertex-star electric field, reconstructed nodal SG carrier currents,
carrier-specific `Eparallel`, and nodal ionization coefficients.

The audit found that the former continuous-cell path diagnostic reconstructed
those quantities a second time. In particular, it replaced the SG-current
direction by a quasi-Fermi-gradient direction when evaluating `Eparallel`.
Consequently the former 7.277244958 V rank-3 crossing did not use the closed
nodal physics.

## Implementation changes

- SG avalanche edge records now retain exact endpoint electric fields,
  electron/hole current vectors, `Eparallel`, and alpha values.
- `continuous_cell` tracing consumes that endpoint ledger directly when
  `source_mapping_mode=nodal_eparallel_p1`; the legacy reconstruction remains
  available for other source modes.
- SG current directions below `1e-7` of the maximum nodal current fall back to
  the electric-field tangent during this benchmark audit.
- Nodal local maxima separated by at most two mesh edges are grouped as one
  physical seed corridor. This prevents the adjacent 327/1381 maxima from
  consuming two `BreakAtIonIntegral` ranks.
- The validation runner now exposes the Eparallel field-recovery and source-
  mapping options needed to reproduce the audit.

## Regression verification

- impact-ionization/node-ledger tests: 51 cases, 725 assertions passed;
- path integration, continuous tracing, stop field, and two-ring seed
  clustering: 9 cases, 41 assertions passed;
- DC sweep and diagnostic CSV tests: 80 cases, 2982 assertions passed.

## Rank-3 crossing

The third physical path remains the 327/1381 drain-junction corridor through
the threshold interval:

| Bias (V) | Representative seed | Rank-3 mean integral |
|---:|---:|---:|
| 8.000 | 1381 | 0.93791856 |
| 8.100 | 1381 | 0.96457611 |
| 8.200 | 1381 | 0.99287531 |
| 8.223 | 1381 | 0.99959954 |
| 8.224 | 1381 | 0.99989383 |
| 8.225 | 1381 | 1.00018828 |
| 8.300 | 1381 | 1.02275456 |

Linear interpolation of 8.224 and 8.225 V gives
`8.224360569 V`. This replaces the obsolete `7.277244958 V` result from the
quasi-Fermi-direction path implementation.

The dense Sentaurus `WriteAll` run has two distinct high-integral corridors at
9.151 V and a third by 9.251 V. Relative to that topology-appearance bracket,
Vela creates and crosses the third path about 1.03 V early. The official
10.448267 V Sentaurus stop remains an adaptive-step overshoot and is not a
dense threshold interpolation.

## Top-three endpoint comparison at 10.4482667308 V

The minimum-distance assignment now maps Vela ranks 1/2/3 directly to
Sentaurus ranks 1/2/3.

| Rank | Vela seed | Vela e | Vela h | Vela mean | Sentaurus mean | Mean ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 990 | 1.515512 | 2.060743 | 1.788127 | 1.794564 | 0.9964 |
| 2 | 1460 | 1.452109 | 1.941337 | 1.696723 | 1.544187 | 1.0988 |
| 3 | 327/1381 | 1.177817 | 1.314998 | 1.246408 | 1.428772 | 0.8724 |

The symmetric mean path-support distances are 0.0819, 0.0208, and 0.0340 um.
Peak electric-field ratios are 1.0188, 0.7498, and 1.0147. On rank 3 the peak
electron/hole alpha ratios are already 1.0133 and 1.0315, while the integral
mean remains 12.8% low. The residual is therefore path length/distributed
support and quadrature geometry, not the local peak ionization model.

## Remaining discrepancy

One global stop-field value cannot close all three paths. Rank 1 currently
extends to the substrate boundary (`y=1 um`) while the imported Sentaurus
plateau support ends near 0.211 um. Conversely, Vela rank 2 and rank 3 supports
are shorter than the corresponding Sentaurus plateau supports. Raising the
stop field would shorten rank 1 but worsen ranks 2 and 3.

This was the next implementation target at the time of this audit.  The later
WriteAll controls established that the numeric path ordering uses the
arithmetic mean of the logged electron and hole injection integrals.  The TDR
`MeanIonIntegral` plateau is an approximate geometry field and must not be used
as the exact BreakAtIonIntegral scalar.  Carrier/path-specific termination and
separate electron/hole tracing remain independent of that plotting field.

## Artifacts

- selected high-voltage branch:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_nodal_colocated_path_final_audit_20260805`;
- dense threshold states:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_nodal_colocated_path_dense_crossing_20260805`;
- endpoint geometry/physics comparison:
  `sg_current_nodal_colocated_path_final_audit_20260805/path_compare/summary.json`.
