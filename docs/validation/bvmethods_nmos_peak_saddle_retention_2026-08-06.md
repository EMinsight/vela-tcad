# BVmethods NMOS peak, saddle, and path-retention audit

## Outcome

The voltage-resolved comparison rejects a fixed electric-field prominence
threshold as the mechanism that turns the p-surface path on. Sentaurus path 3
changes from a weak path to a strong path between 9.151 V and 9.251 V while
its maximum field decreases smoothly from 37.526 MV/m to 37.325 MV/m. Vela's
matching surface maximum behaves the same way: its peak field, merge-tree
saddle, and prominence all vary smoothly through the transition.

The implemented Sentaurus-aligned surrogate therefore separates three
concepts that were previously conflated:

1. all numbered local maxima remain in the WriteAll-style diagnostic output;
2. maxima separated by at most two mesh edges share one physical-corridor ID,
   and the break ranking counts that corridor only once;
3. a boundary minority-carrier quasi-Fermi direction is retained only when its
   seed magnitude is observable relative to the global maximum. The fixed
   relative floor is 5.1e-3; below it, the path uses the electric-field
   direction. There is no bias threshold or seed-node special case.

The exact proprietary Sentaurus peak-retention implementation is not exposed
by the 2018 documentation, so item 3 is an evidence-backed compatibility rule,
not a claim of source-level identity.

## Peak and saddle evidence

The Vela saddle is the highest bottleneck field on a mesh-graph path from the
local maximum to a stronger maximum. Prominence is peak minus saddle.

| bias (V) | seed | peak (MV/m) | saddle (MV/m) | prominence ratio | minority QF relative magnitude | retained mean |
|---:|---:|---:|---:|---:|---:|---:|
| 7.0 | 1462 | 39.039 | 25.645 | 0.3431 | 0.004359 | 1.07e-6 |
| 9.0 | 1462 | 36.468 | 24.335 | 0.3327 | 0.004970 | 5.02e-7 |
| 9.2 | 1462 | 36.128 | 24.344 | 0.3262 | 0.005058 | 3.82e-7 |
| 9.3 | 1461 | 35.960 | 24.293 | 0.3245 | 0.005104 | 1.359452 |
| 10.448267 | 1460 | 33.852 | 22.383 | 0.3388 | 0.005774 | 1.542601 |

The prominence ratio does not rise at the transition. The quasi-Fermi
observability measure crosses the fixed 5.1e-3 floor between 9.2 V and 9.3 V,
which selects the strong minority-carrier trajectory.

Sentaurus shows the corresponding discontinuity without a field discontinuity:

| bias (V) | path | max field (MV/m) | electron | hole | arithmetic mean |
|---:|---:|---:|---:|---:|---:|
| 7.151 | 3 | 40.430 | 1.14e-20 | 0.003429 | 0.0017145 |
| 9.151 | 3 | 37.526 | 7.75e-11 | 0.0018275 | 0.00091375 |
| 9.251 | 3 | 37.325 | 1.28622 | 1.49227 | 1.389245 |

Thus the observed event is a retained-path/channel change, not the birth of a
new total-field maximum.

## Full frozen-branch replay

All 16 frozen states from 7.0 V through 10.4482667308 V replayed successfully.
After collapsing numbered aliases by physical-corridor ID, the relevant ranks
are:

| simulator/bias | physical rank 1 | physical rank 2 | physical rank 3 |
|---|---:|---:|---:|
| Sentaurus 7.151 V | 1.182485 | 1.004786 | 0.0017145 |
| Vela 7.2 V | 1.215399 | 1.035064 | 1.73e-5 |
| Sentaurus 9.151 V | 1.557175 | 1.301600 | 0.00091375 |
| Vela 9.2 V | 1.563784 | 1.213168 | 0.00010317 |
| Sentaurus 9.251 V | 1.575700 | 1.389245 | 1.312110 |
| Vela 9.3 V | 1.581952 | 1.359452 | 1.232268 |

The Vela transition bracket is 9.2--9.3 V and the Sentaurus bracket is
9.151--9.251 V. Their overlap is 9.2--9.251 V. Because the retained path
switches discontinuously, linear interpolation across this bracket is not
physically meaningful.

At the common terminal state, the three distinct physical means are
1.788127, 1.542601, and 1.477063 in Vela versus 1.808320, 1.578535, and
1.468030 in Sentaurus. Relative errors are -1.12%, -2.28%, and +0.62%.

## BreakAtIonIntegral scope

This path threshold is used to stop the non-self-consistent ionization-path
sweep after enough high-field states have been sampled. It is not the
6.377494 V IIC current-intersection result. The sparse official run reaches
10.448 V because its accepted voltage step overshoots the dense path-retention
transition; the current-intersection BV remains a separate comparison target.

## Reproducible artifacts

- audit script: `scripts/audit_bvmethods_nmos_peak_saddle_retention.py`;
- all Vela peaks and saddles: `vela_peak_saddle_by_bias.csv`;
- Vela physical-corridor ranking: `vela_distinct_physical_ranked_by_bias.csv`;
- parsed Sentaurus WriteAll inventory: `sentaurus_writeall_final_by_bias.csv`;
- Sentaurus distinct ranking: `sentaurus_distinct_ranked_by_bias.csv`.

The generated CSV files are under
`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/peak_saddle_retention_audit_20260806`.

## Verification

- `test_path_ionization_integral`: 12 test cases, 78 assertions passed;
- path-ionization subset of `test_dc_sweep`: 3 test cases, 55 assertions passed;
- all 16 frozen postprocess states converged;
- audit script compiles with Python `py_compile` and regenerates the tables.
