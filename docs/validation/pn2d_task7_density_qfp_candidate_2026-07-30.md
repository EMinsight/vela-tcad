# PN2D Task 7 density/QFP support candidate qualification

Date: 2026-07-30

Outcome: `no_authorized_candidate`.

## Candidate and scope

The controlled candidate changes one axis only:

`solver.impact_ionization.quasi_fermi_carrier_truncation = 1.0e-2`

The default remains `0`. When enabled, only the quasi-Fermi potential used by
the avalanche driving field is rebuilt with
`n_eff=max(n,0.01*ni)` and `p_eff=max(p,0.01*ni)`. Physical carrier state,
continuity transport, contacts, mobility drive, source geometry, Van
Overstraeten coefficients, and all production defaults are unchanged.

The `triangle_gss_gradqf_truncated` capability guard now permits this existing
default-off diagnostic so it can be tested without changing a second axis.
Real/AD source-Jacobian coverage includes both `0` and `1.0e-2`.

## Exact-lattice and determinism result

The candidate completed avalanche-off, IIC/postprocess-only, and avalanche-on
at all 29 exact points from `0` through `-20 V`. The process manifest contains:

- 68,034 field records;
- 783 aggregate records; and
- 87 exact-target Newton-attempt records.

Manifest SHA-256:

`94e22044f459791ba43f80931a5ab87293834fe05613fbe912c308826197af1c`

Two independent candidate avalanche-on runs produced identical IV SHA-256:

`4b3bd119e201d60b523b97d5800ef92949fd5706790b30f9d015540bc34dea42`

This is also the baseline IV SHA-256. All curve rows are byte-identical.

## Curve and internal scorecard

| Metric | Baseline | Candidate | Improvement |
|---|---:|---:|---:|
| knee log-current RMSE (dex) | `11.4007360393` | `11.4007360393` | `0%` |
| V_break error (V) | `0.232` | `0.232` | `0%` |
| V_slope error | unavailable | unavailable | none |
| QFP RMSE (V) | `0.3804415239` | `0.3804415239` | `0%` |
| density log-RMSE (dex) | `5.7568129571` | `5.7568129571` | `0%` |
| maximum global error worsening (dex) | n/a | `0` | passes guard |

The baseline and candidate Vela chain inputs each contain 71,733 process
records, and their `records` arrays are exactly equal. The candidate
fixed-transition campaign contains 18 cases and contributes 2,916 Newton
residual/update records.

The WP7 rerun:

- remains `density_qfp_feedback_cause`;
- retains `state` as the first departure at `-19.7/-19.8 V`;
- has no missing stage;
- passes all 203 source/terminal closure rows; and
- does not improve the named internal causal metric.

## Decision gate

Passed guards:

- complete exact lattice;
- duplicate determinism;
- Tasks 4-6/WP7 closure preservation;
- no global-lattice regression above `0.02 dex`; and
- production default unchanged.

Failed authorization gates:

- knee-window RMSE improvement is below 50%;
- neither `V_break` nor `V_slope` improves;
- no nonmonotonic interval is removed; and
- neither QFP nor density internal metric improves.

The machine-readable scorecard is:

`build-release/pn2d-task7-qftrunc1e2-scorecard-20260730/acceptance.json`

Its SHA-256 is:

`29bf422cffc11e1a2d5ebf60fd7cbb2a789fda5427638892ea460fca22cbada4`

Task 8 is not authorized. Further work requires new cross-bias one-stage
density/QFP substitution evidence; increasing the truncation value would be an
empirical scale search and is outside this gate.

## Verification

- focused impact/cell/element/fixed-state/Newton/DC tests: passed;
- 40 focused Python observability tests: passed;
- complete Release CTest: 503/503 passed;
- `git diff --check`: passed.
