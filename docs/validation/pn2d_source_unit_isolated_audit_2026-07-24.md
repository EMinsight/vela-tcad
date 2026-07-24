# PN2D continuity source-unit isolated audit

## Decision

Retain `continuitySourceIntegralFactor() = 1e-8` for TCAD-internal units.
The earlier accurate forward-IV result does not contradict this decision:
that smoke deck is numerically insensitive to the source factor.

## Isolated comparison

Two Release runners were built from the same current source tree.  The only
difference was

`continuitySourceIntegralFactor() = 1`

versus

`continuitySourceIntegralFactor() = (1e-6 m)^2 / (1e-4 m^2/V/s) = 1e-8`.

They were evaluated on the identical imported Minimal6 states at
mirror/sketch x -1/-10/-20 V.

| Check | Result |
|---|---:|
| maximum SG edge-flux relative change | `0.0` |
| maximum SG carrier-term relative change | `0.0` |
| maximum SRH ratio error from `1e-8` | `3.308722e-24` |
| maximum impact ratio error from `1e-8` | `3.308722e-24` |
| analytic/FD Jacobian maximum, factor 1 | `5.048694e-8` |
| analytic/FD Jacobian maximum, factor `1e-8` | `2.989497e-9` |
| first-update pairs reduced by factor-scaled branch | `48/48` |

The denominator `C0*D0` is common to the SG term and the integrated source
term after they are expressed in the same row units.  Dividing both terms by
that denominator cannot cancel a missing relative area/mobility conversion.

## Forward-IV control

The IGBT high-injection smoke sweep converges at all seven points before and
after the patch.  The maximum unit-scaling current change is only
`8.591017e-11` relative.  It is therefore an insensitive regression control,
not evidence that the factor cancels.

The historical legacy-SI and unit-scaling IGBT decks differ by about
`2.8 dex` even before the patch and are not accepted as a strict cross-deck
parity reference.

## Determinism

The isolated `20260724-a` and `20260724-b` roots are byte-identical for every
stable CSV, report, and independent verification JSON.  Both independent
verifications pass.

