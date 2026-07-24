# PN2D Minimal6 Phase D mobility and support alignment

Date: 2026-07-24

Status: `valid`

Primary outcome: `proprietary_model_difference`

Secondary outcome: `support_mismatch`

## Scope

Phase D compared the active Sentaurus and Vela mobility configurations without
changing `include/` or `src/` formulas and without fitting any parameter. The
comparison covers:

- 40 exact states: `mirror/sketch x -1..-20 V`;
- 160 native Sentaurus elements and 320 carrier-element samples;
- 720 carrier-global-edge rows, including 400 rows with an available Vela
  production mobility, 160 geometric-zero rows, and 160 boundary rows without
  a production mobility;
- 960 carrier-element-edge adjacency rows; and
- 160 adjacent-element rows for the central edge `1-5`.

The field-derived region-cell permutation passed with a maximum relative
electric-field residual of `2.959312421212771e-16`.

## D1 - sealed model and parameter contract

The Sentaurus command activates:

- `Mobility(DopingDependence HighFieldSaturation)`;
- `DopingDependence` formula 1, the Masetti form; and
- native element output for electron/hole QFP gradients and mobility.

The Vela deck activates:

- `model: masetti_field`; and
- `high_field_driving_force: quasi_fermi_gradient`.

After converting Sentaurus `cm2/(V s)`, `cm-3`, and `cm/s` parameters to SI,
all 22 numeric carrier parameters match the Vela values. This includes
`mu_const`, `mumin1`, `mumin2`, `mu1`, `Pc`, `Cr`, `Cs`, the Masetti
exponents, the high-field exponents, and the saturation velocities. Both
branches use 300 K. Sentaurus does not expose the internal element
interpolation used to evaluate the high-field mobility.

Therefore the documented Sentaurus parameter substitution is a numeric no-op;
`parameter_mismatch` is false.

## D2 - native element comparison

This replay evaluates the Vela law using the native Sentaurus element
`eGradQuasiFermi` or `hGradQuasiFermi` magnitude, not a Vela triangle
gradient.

| Carrier | Doping control | N | Median absolute error (dex) | P95 (dex) | Maximum (dex) |
|---|---|---:|---:|---:|---:|
| electron | cell-average net doping | 160 | 0.593706 | 0.703886 | 0.724548 |
| electron | arithmetic node mobility | 160 | 0.630953 | 0.733994 | 0.752567 |
| hole | cell-average net doping | 160 | 0.047364 | 0.127339 | 0.184062 |
| hole | arithmetic node mobility | 160 | 0.055805 | 0.147906 | 0.227332 |

The Phase D parity-candidate target is median at most `0.03 dex` and P95 at
most `0.10 dex`. Neither carrier passes. The electron result is especially
diagnostic: the Sentaurus native element electron mobility cannot be replayed
from the exported native element electron QFP-gradient magnitude with the
documented parameters.

The earlier approximately `0.052688 dex` electron element statistic is not a
same-native-field result: it used the Vela triangle QFP-gradient magnitude.
It remains useful as a Vela-support control but is superseded for the Phase D
native-element gate.

## D3 - coefficient-weighted box-edge mobility

For each edge, the Sentaurus-equivalent mobility is

`sum(kappa_element_edge * mobility_element) / sum(kappa_element_edge)`.

The signed dex difference is split exactly into:

1. `native_cell_model`: Vela law on the native Sentaurus element field versus
   Sentaurus native element mobility;
2. `global_support`: Vela production global-edge mobility versus the
   coefficient-weighted Vela native-cell branch; and
3. `production_total`: Vela production global edge versus the
   coefficient-weighted Sentaurus element mobility.

| Carrier | Branch | N | Median absolute error (dex) | P95 (dex) | Maximum (dex) | Absolute-current-weighted mean (dex) |
|---|---|---:|---:|---:|---:|---:|
| electron | native cell model | 200 | 0.599329 | 0.707583 | 0.724548 | 0.572764 |
| electron | global support | 200 | 0.568197 | 0.781317 | 0.819109 | 0.623583 |
| electron | production total | 200 | 0.062732 | 0.961922 | 1.063451 | 0.061720 |
| hole | native cell model | 200 | 0.050123 | 0.133481 | 0.184062 | 0.074522 |
| hole | global support | 200 | 0.139684 | 0.571615 | 0.675545 | 0.139695 |
| hole | production total | 200 | 0.057630 | 0.616202 | 0.710167 | 0.066569 |

The maximum signed-log decomposition closure is
`1.3877787807814457e-16 dex`.

For electrons, the small `0.062732 dex` production-total median is a
cancellation between a large native-element model/interpolation residual and
a large, oppositely signed global-support residual. It must not be interpreted
as native-element mobility parity.

## D4 - central edge 1-5 tail

At -20 V, both topologies give the same mapped result.

| Carrier | Sentaurus adjacent element mobilities (m2/V/s) | Vela native-cell branch (m2/V/s) | Vela local edge (m2/V/s) | Vela global edge (m2/V/s) | Reference current (A/um) | Vela-global candidate (A/um) | Error (dex) |
|---|---:|---:|---:|---:|---:|---:|---:|
| electron | 0.00913328, 0.00913328 | 0.0465705, 0.0484368 | 0.141653 | 0.105701 | 1.625904e-20 | 1.881681e-19 | 1.063451 |
| hole | 0.00684108, 0.00684108 | 0.00763344, 0.00718412 | 0.0470477 | 0.0350988 | 1.546089e-20 | 7.932346e-20 | 0.710167 |

The tail contains both effects:

- native element mobility differs from the Vela law evaluated on the exported
  native element field; and
- one endpoint/global-edge QFP field and mobility replace two element
  mobilities in the Vela production edge operator.

The relative error is large, but the reference carrier-edge current is only
about `1.5-1.6e-20 A/um`.

## D5 - no-fit controls

- Substituting all documented Sentaurus numeric parameters changes the current
  Vela element candidate by exactly zero because the tables already match.
- Algebraically inverting the documented high-field law to obtain an effective
  Sentaurus low-field mobility, then replaying the same high-field law, closes
  all 320 native values to a maximum relative error of
  `3.7683078942700315e-16`.
- The inverted low-field value is diagnostic-only. It is not a fitted
  production candidate and does not identify Sentaurus's internal
  interpolation or hidden driving-force evaluation.

## Decision

Phase D passes its typed exit gate:

- primary: `proprietary_model_difference`;
- secondary: `support_mismatch`;
- `parameter_mismatch: false`; and
- native-element parity target: failed with bounded residuals.

No production mobility formula is justified for modification by this result.
Phase E should carry both a Sentaurus-equivalent coefficient-weighted box-edge
mobility branch and the unchanged Vela production mobility branch into the
fixed imported-state continuity residual.

## Evidence and verification

- Evidence root A:
  `build-release/pn2d-minimal6-phase-d-mobility-support-20260724-a`
- Evidence root B:
  `build-release/pn2d-minimal6-phase-d-mobility-support-20260724-b`
- Independent verifier: passed on both roots.
- A/B directory diff: byte-identical.
- Native inferred-mobility replay maximum relative error:
  `3.7683078942700315e-16`.
- Edge decomposition maximum absolute closure:
  `1.3877787807814457e-16 dex`.
