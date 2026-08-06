# BVmethods NMOS path retention and MeanIonIntegral closure

## Result

The local-peak path topology is now represented by nodal local maxima plus
numbered peak groups.  Nearby aliases retain their Sentaurus-style seed
numbers but share the strongest trajectory in the same two-hop peak group.
At 10.4482667308 V this produces 11 Vela paths, four strong paths above one,
and an exact rank-3/rank-4 duplicate:

| rank | seed node | electron | hole | arithmetic mean |
|---:|---:|---:|---:|---:|
| 1 | 990 | 1.515512276 | 2.060742713 | 1.788127494 |
| 2 | 1460 | 1.452108827 | 1.941337419 | 1.696723123 |
| 3 | 327 | 1.177817387 | 1.314998273 | 1.246407830 |
| 4 | 1381 | 1.177817387 | 1.314998273 | 1.246407830 |

The Sentaurus WriteAll inventory contains ten numbered paths and four strong
paths.  Paths 0 and 1 are identical local-peak aliases.  Sorting its distinct
carrier pairs by `(electron + hole) / 2` gives:

| rank | WriteAll path | multiplicity | electron | hole | arithmetic mean |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 1 | 1.71545 | 1.90119 | 1.80832 |
| 2 | 3 | 1 | 1.39435 | 1.76272 | 1.578535 |
| 3 | 0 | 2 | 1.28569 | 1.65037 | 1.46803 |

The corresponding Vela/Sentaurus mean ratios are 0.9888336, 1.0748720, and
0.8490343.  Rank 1 is therefore within 1.12%.  Rank 2 remains 7.49% high and
rank 3 remains 15.10% low.

## MeanIonIntegral interpretation

Fixed-state Sentaurus controls were run for baseline, electron-only,
hole-only, half-electron, and half-hole avalanche coefficients under both
Eparallel and ElectricField driving forces.  They reject treating the plotted
TDR `MeanIonIntegral` plateau as either the pointwise arithmetic mean of the
two plotted carrier fields or a line integral with simply halved coefficients.

The WriteAll log provides the exact carrier injection integrals used for path
reporting.  Its path order is exactly reproduced by the arithmetic carrier
mean.  The TDR `MeanIonIntegral` field is therefore retained for plateau
support and geometry only; numeric comparison and BreakAtIonIntegral ranking
use the WriteAll carrier values.

## High-voltage branch

All 16 replayed states from 7.0 V through 10.4482667308 V converged.  With the
new numbered peak groups, rank 3 is 0.9928753 at 8.2 V and 1.0227546 at 8.3 V.
This bracket is a Vela diagnostic under the current carrier-integral model; it
is not yet the final Sentaurus voltage closure because rank-3 magnitude is
still 15.10% low at the common final state.

## Implementation and outputs

- Core modes: `nodal_local_maxima` and `numbered_peak_groups` in
  `PathIonizationIntegral`.
- CSV traceability: summary and segment files now export both
  `path_retention` and `seed_mode`.
- Exact-log comparison: `compare_bvmethods_nmos_path_ionization.py` parses the
  final WriteAll inventory, collapses identical carrier pairs with a
  multiplicity count, and uses the arithmetic mean for numeric ranking.
- Corrected comparison output:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/numbered_peak_groups_mean_closed_compare_20260805`.
- Full replay output:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/numbered_peak_groups_nodal_colocated_branch_20260805`.

## Remaining numeric work

Path generation, alias retention, and ordering semantics are structurally
closed.  Remaining work is numeric: shorten or redirect the rank-2 support and
extend the rank-3 carrier-supported trajectory so its electron and especially
hole injection integrals match the WriteAll values.  Avalanche coefficients
and mobility remain fixed during that work.
