# PN2D Minimal6 mobility unit root-cause audit

## Technical summary

The dominant Vela mobility discrepancy is a saturation-velocity unit
conversion defect in `unit_scaling`, not a failure of the Masetti low-field
formula.

The mobility configuration declares saturation velocity in m/s. In the
unit-scaled path, mobility is converted to cm2/(V s) and the QFP-gradient field
is converted to V/cm, but saturation velocity is left numerically in m/s. The
high-field ratio

`low_field_mobility * field / saturation_velocity`

therefore consumes the m/s number as cm/s. This makes the effective physical
saturation velocity 100 times too small and reduces high-field mobility by
approximately 2 dex.

The conclusion is supported by two independent comparisons:

- Across 960 carrier local-edge samples, direct C++ mobility matches the
  unconverted-velocity interpretation to floating-point precision.
- Across 320 native Sentaurus carrier-element samples, correct velocity
  conversion reduces median mobility error from 1.877 to 0.0527 dex for
  electrons and from 1.839 to 0.0478 dex for holes.

No production formula was modified by this audit.

## Exact comparison

| Support | Carrier | Branch | N | Median abs error (dex) | P95 (dex) | Maximum (dex) |
|---|---|---|---:|---:|---:|---:|
| triangle local edge | electron | legacy velocity interpretation | 480 | 9.643e-17 | 1.929e-16 | 2.893e-16 |
| triangle local edge | electron | correct velocity interpretation | 480 | 1.874802 | 1.963048 | 1.968603 |
| triangle local edge | hole | legacy velocity interpretation | 480 | 9.643e-17 | 9.643e-17 | 1.929e-16 |
| triangle local edge | hole | correct velocity interpretation | 480 | 1.789944 | 1.939787 | 1.949354 |
| Sentaurus native element | electron | legacy, cell-average doping | 160 | 1.877046 | 1.947650 | 1.952937 |
| Sentaurus native element | electron | correct, cell-average doping | 160 | 0.052688 | 0.207159 | 0.312971 |
| Sentaurus native element | hole | legacy, cell-average doping | 160 | 1.838943 | 1.923614 | 1.933512 |
| Sentaurus native element | hole | correct, cell-average doping | 160 | 0.047814 | 0.121378 | 0.184062 |

The node-average-doping sensitivity branch gives the same conclusion:

| Carrier | Legacy median (dex) | Correct median (dex) |
|---|---:|---:|
| electron | 1.877004 | 0.058756 |
| hole | 1.838911 | 0.055519 |

## Source-path audit

The production path has four relevant facts:

1. `PhysicalUnitSystem::tcadInternal()` uses micrometers, cm-3, cm2/(V s),
   V/cm, and A/cm2 as internal units.
2. `convertMobilityDefaultsToInternal()` converts Masetti and
   Caughey-Thomas mobility and concentration parameters, but does not convert
   `FieldMobilityParameters::saturationVelocity`.
3. `parseField()` reads JSON keys named
   `electron_saturation_velocity_m_s` and
   `hole_saturation_velocity_m_s` without a `UnitScalingConfig`.
4. `DopingDependentMobility::fieldLimit()` directly evaluates
   `lowFieldMobility * field / params.saturationVelocity`.

Consequently, both the default saturation velocities and explicit JSON
overrides are affected in `unit_scaling`.

## Why support alignment matters

Vela's triangle avalanche path stores mobility on each cell-local edge and
uses the edge QFP difference. Sentaurus native element mobility is a cell
quantity associated with the native element QFP-gradient vector.

Simply averaging Vela local-edge mobility is not the same operator. For a
right-triangle cell, a zero-QFP-difference edge can retain low-field mobility
while the cell QFP-gradient magnitude is high. This can make a physically
correct velocity conversion look worse under an invalid aggregation.

The native-element comparison therefore reconstructs the Vela Masetti/high
field model on the Sentaurus cell field before comparing values. Two doping
interpolation controls were retained:

- arithmetic average of the three nodal mobilities;
- Masetti evaluated at arithmetic cell-average net doping.

Both controls reduce the median error from approximately 1.84-1.88 dex to
approximately 0.048-0.059 dex.

## Method and robustness

- State contract: 40 exact states, mirror/sketch topologies, reverse biases
  -1 through -20 V.
- Direct local-edge contract: 4 cells x 3 local edges x 2 carriers x 40 states
  = 960 rows.
- Native element contract: 4 cells x 2 carriers x 40 states = 320 rows.
- Error metric: absolute base-10 logarithm of the positive mobility ratio.
- The two formula branches differ only by a factor of 100 in saturation
  velocity.
- Independent output roots `20260723-a` and `20260723-b` produced
  byte-identical local-edge, native-element, and summary CSV files.
- An independent CSV-only recomputation reproduced every headline median,
  P95, maximum, and sample count.

## Limitations

- The fixed-state audit directly exposes triangle local-edge mobility, but not
  a separate global SG edge-mobility column.
- Sentaurus native element mobility is not a directed-edge current or flux.
- Residual native-element errors can contain doping interpolation,
  temperature dependence, and solver-specific element evaluation.
- The portable HTML report passed artifact validation and structural
  verification. Browser interaction QA was not run because no installed
  Chromium headless-shell was available; the report builder did not download
  a browser.

## Recommended next steps

1. Add an explicit velocity conversion to the unit system.
2. Apply it to both default and JSON-provided saturation velocities.
3. Add legacy-SI versus unit-scaled mobility parity tests at identical
   physical doping and QFP-gradient fields.
4. Re-run the 40-state fixed-state and self-consistent replacement audits.
5. Resume directed-edge SG current inversion only after mobility parity is
   restored.

The residual question after unit parity is restored is whether the remaining
approximately 0.05 dex median element gap is dominated by doping
interpolation or Sentaurus-specific high-field temperature dependence.
