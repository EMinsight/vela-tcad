# PN2D Task 11 independent scientific review

Date: 2026-07-31

Reviewed snapshot:
`ec003470cb100cea8285b9b05c3396503b0995c8`

Review mode: independent, read-only scientific evidence review. The reviewer
did not modify code, plans, reports, configurations, or generated artifacts
and did not use the independent code-review conclusion.

Verdict:

```text
APPROVE_WITH_CONDITIONS
```

The evidence is sufficient to authorize a separate production-default
proposal and prospective acceptance contract. It is not sufficient to approve
or apply a production-default patch directly.

## Findings

### S1 - balanced M0 raw composite outcome is not pass

Severity: medium; blocks direct default switching.

The machine-readable result at
`build-release/pn2d-task10-balanced-m0-parity-20260731/acceptance.json`
contains:

```text
outcome = ill_conditioned_knee_metric
```

Verified contributing facts:

- Vela `V_break=-19.622 V`;
- Vela `V_slope=-19.828065 V`;
- `abs(V_slope-V_break)=0.206065 V`, slightly above the `0.20 V`
  internal-consistency limit;
- the parity CSV lacks the four closure columns required by the composite
  wrapper, so its closure gate is false;
- the all-domain global maximum is dominated by the `0 V` numerical-current
  floor;
- the all-domain monotonicity gate is false because of the already classified
  low-current state-precision intervals;
- maximum gain error is `0.130837 dex`, above the `0.10 dex` gate.

These failures do not reject the BV-effective high-current agreement, but the
raw composite file must not be described as an unconditional parity pass.

### S2 - M2 lacks a unified machine-readable acceptance artifact

Severity: medium; blocks direct default switching.

The sealed M2 curves, state manifest, determinism evidence, and Task 10 report
support:

- Vela/Sentaurus `V_break=-19.377/-19.391 V`;
- cross-simulator `V_break` difference `0.014 V`;
- median/maximum knee-window current error `0.05663/0.07953 dex`;
- duplicate Vela IV hashes and all `29/29` state hashes match.

Unlike balanced M0, M2 does not have a separately saved, unified
`acceptance.json` produced by the same prospective dual-domain contract.

### S3 - cross-mesh physical convergence is not demonstrated

Severity: medium; mandatory claim boundary.

For M1 to M2:

- Vela `V_break` changes by `0.123 V`;
- Sentaurus `V_break` changes by `0.191 V`;
- both exceed the `0.10 V` mesh threshold;
- maximum integrated-source changes exceed the `2%` threshold by orders of
  magnitude.

The supported claim is same-grid, same-input Sentaurus-golden agreement on
stable M0/M2 comparisons. The evidence does not establish a mesh-independent
physical breakdown knee.

### S4 - M1 is nonblocking only for the same-grid golden objective

Severity: low.

Both simulators exhibit a nonmonotonic M1 avalanche-on branch, so the anomaly
is not a Vela-only same-grid discrepancy. Their jumps occur at different
biases, however, so the evidence does not establish that both solvers follow
the identical nonlinear branch. M1 cannot support a general claim of
coarse-grid robustness or mesh convergence.

### S5 - source, derivative, and deterministic evidence supports the candidate

Severity: none; no blocking finding.

The reviewer confirmed:

- fixed-state integrated-source ratios `1.009537/1.009483`;
- matching-current median and P95 errors below `0.0049 dex`;
- active-node source maximum error below `0.0052 dex`;
- nonzero current-vector direction agreement `100%`;
- branch-resolved Jacobian versus 50-digit finite-difference relative
  difference approximately `1.2e-15`;
- self-consistent QFP RMSE `9.986e-5 V`;
- density log-RMSE `0.001673 dex`;
- Task 7 process-chain closure `203/203`.

The evidence supports the discrete SG/Laux current/source reconstruction and
does not identify a remaining derivative or closure contradiction.

### S6 - no empirical fit or result relocation was found

Severity: none; no blocking finding.

The candidate changes the discrete current/source support through:

```text
current_approximation = element_edge_sg_gss_laux
source_mapping_mode = element_vertex_box_measure
```

It does not alter the Van Overstraeten coefficients, physical parameters,
current scale, voltage axis, Newton tolerances, or empirical field scale. The
midpoint-only negative control reproduces only about `48.6%-48.7%` of the
Sentaurus source, supporting use of the complete SG vector rather than a
post-hoc scale.

## Knee-method review

- `V_slope` is obtained from the first sustained adjacent log-current slope
  above `1 dex/V`, with interpolation.
- `V_break` is the knot from a continuous two-segment linear fit of log
  current, with approximately `1 mV` knot-search resolution.
- `V_curvature` is the location of the maximum adjacent-slope increment.

These are reasonable comparative metrics between paired simulators. They are
not a unique absolute physical definition of breakdown voltage. Fail-closed
handling of a nonmonotonic curve such as M1 is appropriate.

## Conditions before a default patch can be approved

1. Score balanced M0 and M2 with the same prospective dual-domain contract.
2. Bind closure evidence to the exact curve, configuration, mesh, doping, and
   state hashes used by each acceptance result.
3. Save an independent M2 `acceptance.json`.
4. Keep the report boundary explicit:
   same-grid golden agreement is supported; cross-grid physical convergence
   is not demonstrated.
5. Scope the production proposal by device class, mesh quality, compatibility,
   rollback switch, and known M1 topology/continuation risk.
6. Review the resulting default patch separately; this scientific verdict is
   not patch approval.

## Evidence reviewed

- `docs/validation/pn2d_task11_regression_review_2026-07-31.md`
- `docs/validation/pn2d_task10_balanced_mesh_independence_2026-07-31.md`
- `docs/validation/contracts/pn2d_bv_dual_domain_acceptance_v1.json`
- `docs/validation/pn2d_wp6_branch_resolved_jacobian_2026-07-30.md`
- `docs/validation/pn2d_task7_frozen_sg_candidate_2026-07-30.md`
- `build-release/pn2d-task10-balanced-m0-parity-20260731/acceptance.json`
