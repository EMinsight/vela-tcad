# BV performance hotspot localization and improvement

Date: 2026-08-09  
Branch: `codex/bv-performance-profiling`  
Baseline main: `c2f85f6`

## Controlled evidence

The three minimal scenarios use the same physics configuration hash,
`dc28b33f2a95e9664c8d93c7c54b3cb4cac22941175b43f8b36d955a25a2a7a2`:

- high-field transition, `6.08709 -> 6.09959 V`;
- Voltage-to-Current final target;
- external resistor, `1206 V`.

The controlled Release baseline is recorded in
`bv_performance_baseline_c2f85f6_2026-08-08.json`. The independent
RelWithDebInfo `-pg` build and all raw `gmon.out`, flat profiles, call graphs,
and hotspot CSV files are under `build-bv-performance/gprof-971056e`; the
tracked summary is `bv_performance_gprof_971056e_2026-08-08.json`.

MinGW's instrumentation helpers (`_mcount_private` and `__fentry__`) consume
about 42% of sampled time, so absolute profiled wall time is not used for A/B
claims. Stable named self-time hotspots across all three cases are the active
branch fingerprint (5.62-5.82%), Masetti mobility (about 4.25%), the line
search residual callback (2.06-3.03%), the BV process probe (1.75-1.86%), and
SparseLU factorization (1.49-1.70%). MinGW call-graph counts overflow at high
frequency; exact calls come from the internal profiler.

## Internal stage localization

The low-overhead profiler reports stage time and calls plus Newton updates,
residual/Jacobian calls, SparseLU analyze/factorize/solve, boundary full-DD
evaluations, process probes, fingerprints, and continuity diagnostics. The
complete baseline is in
`bv_performance_internal_stages_796c731_2026-08-09.json`.

For the external-resistor case, Jacobian assembly used 2543.0 of 3149.4
seconds (80.75%). The baseline nevertheless spent this cost across seven full
DD solves and 1356 Newton updates, so the implementation first reduced solve
and update counts before changing per-update work.

## Independent optimizations

| Commit | Change | External DD | Newton updates | Wall time (s) |
| --- | --- | ---: | ---: | ---: |
| `796c731` | instrumented baseline | 7 | 1356 | 3149.4 |
| `c496e61` | cross-target state prediction | 7 | 880 | 1961.8 |
| `9ee4f9a` | guarded curvature prediction | 4 | 662 | 1472.8 |
| `943bd9a` | reuse accepted continuity closure | 4 | 662 | 1435.4 |

The final candidate reduces external-resistor full-DD evaluations by 42.9%,
Newton updates by 51.2%, residual calls by 50.5%, Jacobian calls by 51.2%, and
global continuity evaluations from 2735 to 667. Its wall reduction is 54.4%
in the recorded single-run comparison. Because every candidate stage also ran
about 28-31% faster per call than the earlier internal baseline, exact call
reductions are the stable attribution; the close-in-time `9ee4f9a` to
`943bd9a` comparison attributes a further 2.5% wall reduction to continuity
reuse.

## Correctness and remaining limit

The final Voltage-to-Current result is unchanged at `6.395904174606575 V`,
with boundary residual `6.17e-12 A/um`. The external-resistor result is
`6.395887865540998 V`, only `3.18e-11 V` from baseline, with load-line
residual `9.14e-8 V`; global electron and hole continuity ratios are
`3.12e-10` and `3.24e-10`.

The same-current field comparison reproduces the existing validation exactly:
potential P95 error `0.0158109 V`, high-field P95 relative error `17.2860%`,
and high-field correlation `0.992835`. Carrier-density and quasi-Fermi
metrics are also unchanged. The physics hash is unchanged; no mobility,
avalanche, BTBT, or empirical current-scaling parameter was modified.

The external-resistor aspirational limit of at most three full-DD evaluations
is not met: four are required. The guarded predictor reaches a load-line
residual of about `-7.85e-3 V` after the third evaluation, and one final
two-Newton correction is required to reach `9.14e-8 V`. All other stated
performance and correctness gates pass.

The machine-comparable final candidate and the complete three-commit A/B
ladder are in `bv_performance_candidate_943bd9a_2026-08-09.json`.
