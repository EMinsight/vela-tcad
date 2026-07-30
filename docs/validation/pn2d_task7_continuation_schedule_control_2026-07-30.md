# PN2D Task 7 continuation-schedule branch-invariance control

Date: 2026-07-30

Outcome: `continuation_invariant_cross_block_reversal`.

Task 7 outcome: `no_authorized_candidate`.

## Scope

This control changed only the predeclared continuation step schedule. It did
not change the mesh, contacts, doping basis, mobility, recombination, impact
ionization model, avalanche driving force, Van Overstraeten coefficients,
Newton budget, tolerances, source assembly, residual, Jacobian, or production
defaults.

The schedules were:

| Schedule | Initial step | Maximum step | Minimum step | Growth/shrink |
|---|---:|---:|---:|---:|
| `standard_0p05` | 0.05 V | 0.05 V | `1e-10 V` | 1.2 / 0.5 |
| `refined_0p025` | 0.025 V | 0.025 V | `1e-10 V` | 1.2 / 0.5 |

The common physics configuration SHA-256 is
`291f28c1c8e7705301be27a34c532c5cdfa9760945725b1881fbd235b2b29d2c`.
The common configuration hash after removing only schedule and output-path
fields is
`f72b76112a0df481a7da2561ddb2a40496d3f00dd1167a5a71051a7259c26066`.

## Determinism and exact-lattice result

Each schedule was run twice from an independent output root. All four
avalanche-on runs completed all 29 exact target biases through `-20 V`.
Within each schedule, the IV CSV and every exact-state SHA-256 matched between
the duplicate runs.

The two schedules reached the same knee branch within the frozen state gates:

| Bias | max `psi` difference | max QFP difference | max density log difference |
|---|---:|---:|---:|
| -19.7 V | `7.71e-13 V` | `7.96e-13 V` | `1.37e-12 dex` |
| -19.8 V | `1.23e-12 V` | `1.28e-12 V` | `2.84e-12 dex` |

Across all 29 exact targets, the schedule-to-schedule absolute log-current
difference has RMSE `8.63e-4 dex` and maximum `0.0045211 dex` at `-2 V`,
inside the already frozen `0.02 dex` global-worsening limit.

## Named internal causal metric

Both schedules were evaluated with the same Task 6 frozen-state QFP
substitution probe and duplicate-run check.

| Bias | Schedule | carrier-only QFP improvement | full-coupled QFP improvement | carrier/full update cosine |
|---|---|---:|---:|---:|
| -19.7 V | standard | 13.128% | -7.023% | 0.6412 / -0.1994 |
| -19.7 V | refined | 13.128% | -7.023% | 0.6412 / -0.1994 |
| -19.8 V | standard | 12.942% | -7.285% | 0.6367 / -0.2132 |
| -19.8 V | refined | 12.942% | -7.285% | 0.6367 / -0.2132 |

The refined schedule therefore does not remove the carrier-only versus
full-coupled Poisson-QFP update reversal at either adjacent bias. A complete
candidate curve campaign is not authorized.

## Decision

- `task7_outcome`: `no_authorized_candidate`;
- `complete_curve_campaign_authorized`: false;
- `task8_authorized`: false;
- production defaults changed: false;
- next gate: `retain_task8_stop_no_continuation_candidate`.

The remaining observed behavior is schedule-invariant at the tested
resolution. Continuation is not an evidence-authorized correction axis.
Further work requires a separately authorized observation-only decomposition
of the coupled Poisson-QFP cross block; it must not enter Task 8 or change a
production default.

Primary scorecard:

`build-release/pn2d-task7-continuation-schedule-scorecard-20260730/acceptance.json`.

SHA-256:
`d33ac0e6984f0a2f9745e39593cf16b076a53640a727e39838cce2be7f8ab234`.

## Code and tests

- `scripts/run_pn2d_bv_exact_lattice_process.py`;
- `scripts/analyze_pn2d_bv_continuation_schedule_control.py`;
- `tests/regression/test_pn2d_bv_exact_lattice_process.py`;
- `tests/regression/test_pn2d_bv_continuation_schedule_control.py`.

Forty focused process-observability tests and the complete Release CTest suite
(`504/504`) pass.
