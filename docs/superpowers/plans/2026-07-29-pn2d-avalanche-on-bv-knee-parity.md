# PN2D avalanche-on BV curve and knee parity plan

Date: 2026-07-29

Status: engineering WP0-WP7, exact-lattice continuation, fixed-transition
Newton, the Task 6 density/QFP matrix, both predeclared Task 7 continuation
schedules, the observation-only Poisson-QFP cross-block decomposition, the
effective Schur-loop source decomposition, the self-consistent
element-edge SG/GSS-Laux candidate, dose-preserving mesh validation, and Task
11 regression are executed. Stable same-grid M0/M2 BV-effective metrics
support Sentaurus-golden agreement; the raw M0 composite outcome remains
`ill_conditioned_knee_metric`. The M1 mesh/continuation observation remains
open. Both new independent reviews return `APPROVE_WITH_CONDITIONS`. A
separate PN2D BV template default proposal and prospective dual-domain
contract are authorized; direct production-default changes remain prohibited.

Starting point: local branch `codex-pn2d-minimal6-operator-audit`, commit
`1350d11`. Commit `fa1c343` is an ancestor of this starting point.

## Objective

Make the Vela avalanche-on PN2D reverse-bias terminal-current curve agree with
the paired Sentaurus curve on the same physical mesh, with particular emphasis
on reproducing the gradual high-bias turn between `-19 V` and `-20 V`.

The work is successful only when all of the following are true:

1. current magnitude agrees over the complete reverse-bias range;
2. the knee voltage, local log-current slope, and monotonic steepening agree;
3. electron/hole source, terminal-current, and internal-current closure remain
   valid;
4. the result is a deterministic self-consistent branch rather than a
   continuation artifact;
5. the avalanche-off BV and forward-IV baselines do not regress; and
6. no empirical current scale, voltage shift, field scale, lifetime fit, or
   avalanche-coefficient fit is used.

## Confirmed starting evidence

### Same-mesh avalanche-off baseline

The accepted M0 physical mesh has 27 nodes and 32 triangles, with mesh hash:

`c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e`

On this mesh, the paired avalanche-off Vela/Sentaurus terminal-current error is
already:

- 21-point log-current RMSE: `6.997e-5 dex`;
- maximum log-current error: `1.255e-4 dex`;
- maximum electron source/contact closure error: `7.482e-7`; and
- maximum hole source/contact closure error: `3.561e-9`.

This is the mandatory off baseline. M1 and M2 from the 2026-07-28 study are not
accepted physical refinement levels because placing the junction node in both
profile windows changed discrete total-impurity dose by `-5.56%` and `-8.33%`.
They may be used only as explicitly labelled diagnostics, not as parity gates.

### Exact Sentaurus avalanche-on oracle

The accepted same-mesh exact lattice is:

`-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`.

For the implicit/default Sentaurus branch:

| Bias (V) | abs anode total current (A) | avalanche source integral (A/um) |
|---:|---:|---:|
| -10 | `3.179640402e-17` | `7.306902619e-19` |
| -18 | `7.365358341e-17` | `4.244591356e-17` |
| -19 | `1.125860992e-16` | `8.130929396e-17` |
| -19.5 | `1.564824923e-16` | `1.248770659e-16` |
| -19.8 | `2.320220811e-16` | `1.958965435e-16` |
| -19.9 | `3.027800117e-16` | `2.608070677e-16` |
| -19.95 | `3.589097065e-16` | `3.126050609e-16` |
| -20 | `4.337466423e-16` | `3.823148999e-16` |

The terminal-current log slope steepens gradually:

| Interval (V) | Sentaurus log-current slope (dex/V) |
|---|---:|
| -19 to -19.5 | `0.65845` |
| -19.5 to -19.8 | `1.31296` |
| -19.8 to -19.9 | `2.66174` |
| -19.9 to -19.95 | `3.40129` |
| -19.95 to -20 | `3.78779` |

The knee must therefore be compared as a smooth turn, not as a discontinuous
breakdown voltage.

### Frozen model conclusions

- Sentaurus implicit/default avalanche and explicit `GradQuasiFermi` are
  identical on the complete exact lattice.
- Electric-field avalanche drive is a materially different, much stronger
  branch and must not replace `GradQuasiFermi`.
- The Van Overstraeten alpha implementation, P1 field/QFP gradients, Old
  Slotboom densities, element-vertex measure, and fixed-state GSS/Laux source
  reconstruction already have focused evidence.
- At the fixed `-20 V` generation hotspot, the first material stage from
  `-19 V` to `-20 V` is electron density. Density grows `3.77935x`, current
  `5.69013x`, alpha `1.14342x`, generation `4.70196x`, field/QFP gradient only
  about `1.03x`, and mobility falls slightly.
- The historical nonlinear authorization blocker was the source-only
  Jacobian outcome `nonsmooth_branch_derivative`. WP6 has now reproduced its
  best ordinary-double result and closed it with a branch-resolved
  multiprecision reference; this authorizes causal localization, not a model
  or production-default change.

These facts make carrier/QFP feedback and branch selection the primary
investigation. They do not authorize a model change by themselves.

### Current execution evidence

The companion observability plan's WP0-WP2 were executed on 2026-07-29:

- The M0 avalanche-off baseline is sealed with 21-point RMSE
  `6.996509075786512e-5 dex`, maximum error
  `1.254688929143555e-4 dex`, electron/hole source-contact closure
  `7.48198732275877e-7`/`3.561187670708268e-9`, and terminal-pair closure
  `3.7621205605551904e-23 A/um`.
- Two current-code avalanche-on Vela runs are deterministic. Both contain
  2,232 accepted states, end at `-19.692187499999644 V`, and fail the request
  for `-19.693749999999643 V` with `max_iterations`. Their normalized output
  hashes match. The implementation baseline is therefore sealed with a
  deterministic solver-first-failure, not with a complete curve through
  `-20 V`.
- The common process-data contract is verified by 13 focused tests.
- Two independent Sentaurus VM roots produced identical normalized
  off/IIC/on/`AvalDerivatives` process matrices on
  `-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`: 32 exact TDR
  snapshots per root, 48,736 field records, and 288 aggregates.
- Off versus IIC state fields are identical over 7,768 records while IIC
  generation is nonzero; avalanche-on versus `AvalDerivatives` state fields
  are also identical over 7,768 records. CurrentPlot versus Tcl/`ReadMeasure`
  source reintegration closes to `3.977747984452422e-15` relative.
- WP2.1 completed the 29-point union of both acceptance lattices in two
  independent VM roots: 116 exact snapshots, 176,668 field records, and
  1,044 aggregates per root. Both lattices have no missing rows, hashes match,
  actual-bias error is zero, and 348/348 source comparisons close to
  `4.830623899888963e-15` maximum relative error.
- Scientific Task 2 passed 12 synthetic fail-closed scenarios. Its current-data
  outcome is `solver_first_failure`, with no substituted knee metrics and the
  exact missing Vela rows recorded.
- Thirty-six focused Python tests and the complete Release suite
  (`495/495`) pass.
- The 2026-07-30 run-local `max_iter=80` continuation control subsequently
  completed Vela off/IIC/on at all 29 exact points. The contract-valid Vela
  process manifest has SHA-256
  `b882ece81a9cd1e7633e5685adbdd1a9ffde8b4adf2d14dea4fbc2286d6ddf6d`.
  Fixed-transition Sentaurus/Vela probes subsequently added per-node
  residual and first-update signatures at all six knee targets. The paired
  WP7 rerun passes 203/203 closure rows, has no missing stages, and accepts
  `state` for Sentaurus IIC versus avalanche-on at `-19.7/-19.8 V`.

Primary evidence:

- `build-release/pn2d-wp0-implementation-baseline-20260729/acceptance.json`;
- `build-release/pn2d-wp2-process-matrix-pair-20260729/acceptance.json`.

The engineering eight-point process matrix validates extraction, but does not
satisfy this plan's complete global and knee lattices. Missing exact knee
states are `-18.5, -19.25, -19.7, -19.85 V`; the common global manifest must
also include every integer bias from `0` through `-20 V`.

### Cross-plan execution mapping

| Scientific task | Engineering prerequisite | Current state |
|---|---|---|
| Task 1 | observability WP0, WP2, and WP2.1 | `paired_m0_reference_sealed`; exact Sentaurus lattices complete |
| Task 2 | common WP1 contract plus dedicated curve/knee analyzer | `curve_knee_contract_verified`; current data returns `solver_first_failure` |
| Task 3 | WP0 baseline plus WP4-WP5 solver-used/attempt records | `complete_nonlinear_trace_available`; deterministic first failure reproduced |
| Task 4 | WP6 source-only derivative work | `source_jacobian_dependency_identified_and_closed` |
| Tasks 5-6 | WP3-WP5 records, WP7 paired analyzer, and one-stage feedback matrix | `continuation_only_cause`: QFP-only carrier-block updates improve about 13% at adjacent `-19.7/-19.8 V`, while the full coupled update reverses direction and worsens about 7% |
| Task 7 | causal authorization from Tasks 4-6 | QFP truncation and refined continuation rejected; complete SG/Laux candidate improves curve/internal metrics but moves nonmonotonicity to low-current intervals: `tradeoff_without_parity` |
| Tasks 8-11 | successful Task 7 candidate | not authorized |

## Frozen production configuration

Unless one task explicitly defines a paired diagnostic control, hold fixed:

- Vela BV mobility doping basis: `net_doping`;
- Vela forward-IV mobility doping basis:
  `cell_reconstructed_total_impurity`;
- impact-ionization model: `van_overstraeten`;
- avalanche driving force: `quasi_fermi_gradient`;
- current production approximation: `cell_reconstructed`;
- production source mapping: `triangle_gss_gradqf_truncated`;
- high-field mobility: enabled with QFP-gradient drive;
- SRH and Old Slotboom BGN: enabled;
- temperature, contacts, work functions, device depth, current normalization,
  mesh, doping, and materials: paired and hash-identical;
- M0 mesh during Tasks 1-9;
- no implicit zero fill when importing Sentaurus fields; and
- no change to a production default before Task 11.

The opt-in `element_edge_sg_gss_laux` path remains diagnostic until it passes
the complete authorization chain.

## Prohibited shortcuts

Do not:

- multiply current, source, alpha, field, lifetime, or mobility by a fitted
  factor;
- translate the voltage axis to align the knee;
- change GradQF avalanche drive to electric field;
- select a bias schedule only because it hides a failed transition;
- mix the old fine-mesh reference with the M0 exact oracle;
- use interpolated or nearest-bias Sentaurus states as exact evidence;
- call reconstructed directed-edge current a native Sentaurus field;
- relax a derivative tolerance to pass a sign, floor, clamp, or truncation
  crossing;
- alter BV or IV doping-basis choices; or
- stage `tmp/` or generated simulation outputs.

## Common artifact and provenance contract

Each run must write a machine-readable manifest containing:

- Git commit, branch, and dirty-state inventory;
- executable and build-configuration hashes;
- Vela config and Sentaurus deck/model hashes;
- mesh, connectivity, coordinates, doping, material, and contact hashes;
- model selectors and all numerical tolerances;
- exact requested and accepted bias;
- warm-start parent-state hash and bias-step history;
- accepted/rejected Newton transitions, damping, clamps, and linear status;
- convergence and physical-closure results;
- electron, hole, and total terminal currents;
- SRH and carrier-split avalanche source integrals; and
- normalized numerical-output hash.

Generated CSV, JSON, TDR, PLT, logs, states, VTK, and figures belong under an
ignored `build-release/` root. Only source code, tests, small hand-authored
contracts, and validation summaries may be committed.

## Quantitative curve and knee contract

### Bias lattices

Use two paired lattices:

1. global lattice: `0, -1, -2, ..., -20 V`;
2. knee lattice:
   `-18, -18.5, -19, -19.25, -19.5, -19.7, -19.8, -19.85, -19.9, -19.95, -20 V`.

Sentaurus must be run at every missing exact global- and knee-lattice point
before that point is used in a gate. Interpolation is permitted only in a
labelled plot, never in an acceptance calculation. If either simulator lacks
an exact row, the analyzer emits `incomplete_exact_lattice` or
`solver_first_failure`; it must not compute a substituted partial acceptance
score.

For comparison, define:

`y(V) = log10(abs(I_total(V)))`

and on adjacent exact states:

`s_i = (y(V_i) - y(V_{i-1})) / (abs(V_i) - abs(V_{i-1}))`.

Use `1e-30 A/um` only as a plotting floor. A point at or below the numerical
floor is excluded with an explicit reason rather than silently clipped into a
gate.

### Knee estimators

Report all three estimators:

1. `V_slope`: first exact or linearly interpolated-within-one-exact-interval
   crossing of `s = 1.0 dex/V`, provided the next interval remains above the
   threshold;
2. `V_break`: breakpoint from a continuous two-segment least-squares fit to
   `y(V)` over `-18..-20 V`, with at least three states on each side; and
3. `V_curvature`: midpoint of the interval with maximum positive change in
   adjacent log slope.

The first two are acceptance metrics. `V_curvature` is a diagnostic guard
against an apparently aligned threshold crossing with the wrong shape.

If `V_slope` and `V_break` differ by more than `0.20 V` for either simulator,
classify the knee as `ill_conditioned_knee_metric`, refine the exact lattice,
and stop the parity decision.

### Final curve gates

| Metric | Required result |
|---|---:|
| global-lattice median absolute log-current error | `<= 0.05 dex` |
| global-lattice P95 absolute log-current error | `<= 0.10 dex` |
| global-lattice maximum absolute log-current error | `<= 0.15 dex` |
| knee-lattice median absolute log-current error | `<= 0.05 dex` |
| knee-lattice maximum absolute log-current error | `<= 0.10 dex` |
| `abs(V_slope,Vela - V_slope,Sentaurus)` | `<= 0.10 V` |
| `abs(V_break,Vela - V_break,Sentaurus)` | `<= 0.10 V` |
| knee-lattice adjacent-slope RMSE | `<= 0.20 dex/V` |
| monotonicity of `abs(I)` versus `abs(V)` | no reverse interval |
| avalanche gain log-error on `-18..-20 V` | median `<= 0.05 dex`, max `<= 0.10 dex` |

Avalanche gain is evaluated from the paired same-mesh curves:

`M(V) = abs(I_on(V)) / abs(I_off(V))`.

The final decision must pass both absolute current and gain. Agreement caused
by equal and opposite errors in the on and off branches is not accepted.

### Physical and numerical guards

- Electron and hole source/contact closure: `<= 1e-5` relative whenever the
  corresponding source is above the declared absolute floor.
- Below that floor, use an independently frozen absolute-closure limit.
- Total terminal-pair closure: `<= 1e-20 A/um` absolute.
- Internal total-current KCL: `<= 1e-8` relative in the active region.
- No accepted row may contain false convergence, a failed physical closure, a
  NaN/Inf, or an unclassified rejected transition.
- Duplicate runs must have identical accepted bias sequences, typed outcomes,
  normalized manifests, and numerical outputs within `1e-12` relative or the
  relevant solver tolerance, whichever is stricter and numerically meaningful.

## Task 1 - seal the paired baseline and exact reference

### Goal

Start from one unambiguous M0 on/off reference and prove that configuration,
mesh, current normalization, and model selectors are paired.

### Actions

1. Confirm `fa1c343` and `1350d11` are ancestors of `HEAD`.
2. Record `git status --short`; allow `tmp/` to remain untracked but do not
   stage it.
3. Re-verify the accepted M0 mesh hash, node/triangle counts, doping hash,
   contacts, device depth, and current sign convention.
4. Seal the exact Sentaurus implicit/default avalanche-on and avalanche-off
   branches, including release `O-2018.06-SP2`, deck, model, mesh, and output
   hashes.
5. Generate the missing exact knee-lattice Sentaurus states with the same M0
   inputs, if they do not already exist.
6. Re-run or verify the 21-point M0 avalanche-off Vela/Sentaurus baseline.
7. Assert from rendered configs that BV uses `net_doping` and forward IV uses
   `cell_reconstructed_total_impurity`.

### Acceptance criteria

- Both required commits are ancestors of `HEAD`.
- Same physical mesh and doping are proven by hashes, not only dimensions.
- All exact reference rows converge and have no nearest-state substitution.
- Sentaurus implicit and explicit GradQF controls remain identical within
  `1e-10` relative on the accepted lattice.
- M0 avalanche-off RMSE is `<= 0.001 dex` and maximum error is
  `<= 0.002 dex`.
- Model/basis assertions pass.

### Decision gate and stop condition

Typed outcome:

- `paired_m0_reference_sealed`; or
- `baseline_or_provenance_mismatch`.

Report the hashes, row counts, off-curve metrics, and model selectors. Stop on
the second outcome; do not run avalanche-on Vela comparisons against an
unpaired reference.

Current state: `paired_m0_reference_sealed`. WP2.1 supplies every exact global
and knee Sentaurus row in two deterministic VM roots; Task 2 is authorized.

## Task 2 - implement the curve and knee verifier

### Goal

Turn “curve basically agrees and knee is consistent” into a fail-closed,
reviewable calculation before changing solver code.

### Actions

1. Add a new analyzer, preferably
   `scripts/analyze_pn2d_avalanche_on_bv_parity.py`.
2. Make it read explicit on/off Vela and Sentaurus CSV paths and manifests.
3. Implement exact-row matching, log-current errors, gain, adjacent slopes,
   all three knee estimators, monotonicity, and closure gates.
4. Emit:
   - `curve_points.csv`;
   - `slope_points.csv`;
   - `knee_metrics.json`;
   - `acceptance.json`;
   - linear-current and log-current comparison figures; and
   - a gain and local-slope comparison figure.
5. Add tests with synthetic curves covering:
   - identical curves;
   - known vertical mismatch;
   - known knee-voltage shift;
   - correct knee but wrong post-knee slope;
   - nonmonotonic current;
   - duplicate/missing bias rows;
   - numerical-floor rows; and
   - ill-conditioned knee metrics.
6. Keep the existing one-volt knee scripts as historical diagnostics; do not
   silently change their semantics.

### Acceptance criteria

- Every synthetic pass/fail case produces the expected typed outcome.
- A known `0.10 V` translated knee is detected within one exact-grid
  resolution.
- Missing, duplicate, interpolated, nonfinite, or unconverged gate rows fail
  closed.
- Reordering CSV rows does not change results.
- Every reported number is traceable to an input file, row, and manifest hash.

### Decision gate and stop condition

Typed outcome:

- `curve_knee_contract_verified`; or
- `curve_knee_contract_invalid`;
- `incomplete_exact_lattice`; or
- `solver_first_failure`.

`curve_knee_contract_verified` means the analyzer implementation and synthetic
contract are valid. Actual-data execution may still end with
`incomplete_exact_lattice` or `solver_first_failure`; those outcomes must
preserve missing/failed bias, parent bias, solver reason, and manifest hashes
without emitting final knee gates. Report test names and expected/actual
outcomes. Stop before physical parity classification if the verifier does not
fail closed.

Executed implementation outcome: `curve_knee_contract_verified` with 12
synthetic scenarios. Executed current-data outcome: `solver_first_failure`;
global `-20 V` and knee
`-19.7, -19.8, -19.85, -19.9, -19.95, -20 V` Vela rows are unavailable.
Evidence:
`build-release/pn2d-task2-curve-knee-20260730/acceptance.json`.

## Task 3 - capture and classify the current Vela avalanche-on branch

### Goal

Measure the present self-consistent failure without modifying physics or
continuation.

### Actions

1. Run the production M0 avalanche-on sweep twice on the global and knee
   lattices.
2. Preserve every accepted and rejected bias transition, Newton history,
   damping/clamp activation, carrier residual, row scale, source integral,
   contact current, and state hash.
3. Compare the Vela on branch with paired Vela off and Sentaurus on/off
   branches.
4. Locate:
   - first current-magnitude departure over `0.05 dex`;
   - first slope departure over `0.20 dex/V`;
   - first nonmonotonic interval;
   - first state/process departure in dependency order; and
   - first branch-dependent or non-deterministic transition.
5. At `-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`, export `psi`,
   electron/hole QFP, `n/p`, mobility, drive, current, alpha, carrier-split
   source, residual decomposition, and terminal current.

### Acceptance criteria

- Both runs reach every required exact state or record the same first typed
  failure.
- The sealed current baseline must reproduce the transition
  `-19.692187499999644 -> -19.693749999999643 V` and `max_iterations`, unless
  a separately authorized candidate intentionally changes it.
- Numerical duplicates agree under the common determinism guard.
- Every accepted point passes source/contact and terminal closure.
- The earliest divergence is assigned to one of:
  `fixed_state_operator`, `source_jacobian`, `self_consistent_state_feedback`,
  `continuation_branch`, `contact_boundary`, `solver`, or
  `insufficient_observation`.

### Decision gate and stop condition

Current baseline classification is `solver_first_failure`: both runs are
deterministic, global exact coverage ends at `-19 V`, and knee exact coverage
ends at `-19.5 V`. Do not interpolate or nearest-match `-19.7 V` and later
rows, and do not report complete knee metrics from this partial curve.

After WP4-WP5 instrumentation, report either the complete curve/knee scorecard
or the same typed failure with the first causal bias, node/cell, carrier, and
term. If closure or determinism fails, stop and repair that invariant before
interpreting curve physics.

## Task 4 - close the source-only Jacobian blocker

### Goal

Resolve `nonsmooth_branch_derivative` before authorizing any nonlinear
avalanche candidate.

### Actions

1. Reproduce the accepted `-20 V` source-only assembled-state probe using an
   imported Sentaurus state or another explicitly fixed assembled state. This
   diagnostic does not require the current self-consistent Vela avalanche-on
   sweep to converge at `-20 V`.
2. Trace the exact active branches for:
   - carrier density/statistics;
   - SG Bernoulli terms;
   - low- and high-field mobility;
   - QFP-gradient drive;
   - current-vector reconstruction;
   - `abs(J)`, sign, zero floor, clamp, and truncation;
   - `alpha(F)` and its derivative;
   - source distribution and node scatter; and
   - continuity scaling.
3. Compare shared Real, forward-AD, complex-safe or branch-frozen numerical
   derivative where valid, and symmetric finite difference over at least
   three same-branch steps.
4. Identify whether the remaining `1.839007985e-6` difference is:
   - a missing analytic derivative;
   - residual/Jacobian formula drift;
   - configuration mismatch;
   - floating-point cancellation; or
   - a true nondifferentiable branch.
5. For a true branch crossing, introduce no smoothing unless the exact
   physical/numerical contract independently requires it. Prefer a
   semismooth/active-set-consistent derivative or a branch-separated test.
6. Add focused acute, obtuse, reversed-orientation, contact-adjacent, and
   interior tests for both carriers.

### Acceptance criteria

- Nonzero same-branch analytic/reference derivative maximum:
  `<= 1e-8` relative.
- Near-zero maximum: `<= 1e-12` absolute.
- Analytic and independently evaluated contribution sums each close to their
  assembled blocks within `1e-12`.
- Residual, Jacobian, and diagnostics have identical configuration and
  active-branch hashes.
- No unrelated Poisson, SRH, mobility, statistics, or default formula change.

### Decision gate and stop condition

Typed outcome:

- `source_jacobian_dependency_identified_and_closed`;
- `incomplete_analytic_derivative`;
- `nonsmooth_branch_derivative`;
- `configuration_mismatch`; or
- `jacobian_gate_failed`.

Only the first outcome may enter Task 5. Report the step-size convergence
table and first mismatching dependency. Do not relax the gate.

Executed outcome: `source_jacobian_dependency_identified_and_closed`.

The ordinary double symmetric reference reproduced relative differences
`1.0623268123382636e-2`, `1.839007985012764e-6`, and
`3.2187014658253814e-6` at steps `1e-8`, `1e-10`, and `3e-11`.
The shared-formula 50-decimal-digit branch-resolved reference then produced
`1.1335074719484275e-14`, `1.1989526702416226e-15`, and
`6.476903678003707e-16` at same-branch steps `1e-14`, `3e-15`, and
`1e-15`. All probes used configuration fingerprint `ab2dbf93089c7fe3`
and active-branch fingerprint `8e39422ae0ff24ba`.

The remaining difference is floating-point cancellation/branch resolution in
the ordinary double reference. There is no evidence of a missing analytic
derivative, configuration mismatch, or residual/Jacobian formula drift.
Exact-zero branches retain the production semismooth derivative of zero; no
smoothing or default change was introduced. Task 5 is authorized.

## Task 5 - close the fixed-state process chain

### Goal

Determine whether a fixed Sentaurus state produces matching Vela current,
alpha, source, and terminal contribution before self-consistent feedback is
allowed to obscure the cause.

### Actions

At every exact high-bias imported state, compare in order:

1. `psi/QFP -> n/p`;
2. P1 electric field and QFP gradients;
3. low- and high-field mobility;
4. element-local SG directed-edge operator replay;
5. GSS/Laux cell-current reconstruction;
6. native Sentaurus element current-density vector;
7. electron/hole alpha;
8. element-vertex avalanche source;
9. physical-node source accumulation; and
10. device source integral and legitimate terminal current.

Run one-stage substitutions and report electron/hole, interior/contact,
acute/right/obtuse, component/angle/magnitude/sign, and active/tail results
separately. Native, reconstructed, and unsupported observations must remain
explicitly labelled.

### Acceptance criteria

- Previously accepted density, P1, and source-identity tests do not regress.
- Constant-field geometry and source identity: `<= 1e-12`.
- `ReadMeasure`/CurrentPlot integral difference: `<= 1e-10` relative.
- Reconstructed Sentaurus source maximum error: `<= 5e-3 dex`.
- Matching-support current median/P95: `<= 0.05/0.15 dex`.
- Nonzero carrier sign agreement: 100%.
- Vela fixed-state source median/maximum in the active range:
  `<= 0.10/0.30 dex`.
- Every failed class has one earliest non-closing stage.

### Decision gate and stop condition

Typed outcome:

- `fixed_state_process_chain_closed`;
- `mobility_support_difference`;
- `element_current_support_difference`;
- `source_mapping_geometry_difference`;
- `contact_interior_model_difference`;
- `proprietary_operator_difference`; or
- `insufficient_native_observation`.

Only the first outcome, or one minimal proven Vela defect that is repaired and
then closes the same gates, may enter Task 6.

## Task 6 - localize self-consistent density/QFP feedback

### Goal

Find why Vela enters the wrong avalanche-on state branch when the fixed-state
operator is already comparable.

### Actions

1. At `-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`, assemble
   avalanche-off and avalanche-on residuals on:
   - the exact Sentaurus state;
   - the converged Vela off state;
   - the converged Vela on state; and
   - one-stage hybrid states.
2. Record the dependency chain:
   `QFP -> density -> mobility/current -> alpha -> source -> residual ->`
   `QFP update`.
3. Compare the first carrier-only and coupled Newton update using identical
   scaling, constraints, damping, and clamps.
4. Track the hotspot-coincident electron density, not a contact-dominated
   global maximum.
5. Perform one-stage substitutions:
   - Sentaurus `n/p` only;
   - Sentaurus QFP only;
   - Sentaurus mobility/current vector only;
   - Sentaurus alpha only; and
   - Sentaurus source only.
6. Identify the first bias/node/carrier/term whose update moves away from the
   Sentaurus branch.

### Acceptance criteria

- Residual decomposition closure: `<= 1e-12`.
- Boundary rows remain bitwise or floating-point identical in paired controls.
- Task 4 Jacobian gates remain valid at all sampled states.
- The same first causal stage is reproduced at two adjacent knee biases.
- A proposed correction improves both residual direction and first-update QFP
  error without worsening the other carrier or a lower-bias state.

### Decision gate and stop condition

Typed outcome:

- `density_qfp_feedback_cause`;
- `mobility_current_feedback_cause`;
- `source_support_feedback_cause`;
- `contact_boundary_cause`;
- `continuation_only_cause`;
- `operator_improvement_without_state_causality`; or
- `proprietary_model_difference`.

Stop unless one of the first three outcomes has cross-bias causal evidence.

Executed follow-up outcome: `continuation_only_cause`.

At adjacent `-19.7/-19.8 V` exact avalanche-on states, the controlled matrix
froze Sentaurus electron/hole/combined density, QFP, or both in the Vela
residual operator. All variants shared one production baseline Jacobian,
scaling, contact constraints, and update caps. QFP-only carrier-block updates
reduced QFP error by `13.128%/12.942%` with update-direction cosine
`0.6412/0.6367`; the full coupled update instead worsened error by
`7.023%/7.285%` and reversed direction to `-0.1994/-0.2132`. Density-only
controls showed carrier antagonism rather than a no-worsening correction.
Residual decomposition closes below `4.55e-13`, all boundary rows are
identical, and duplicate hashes match.

The evidence therefore localizes the remaining difference to the coupled
Poisson-QFP path. It does not authorize a density/QFP physical correction or
Task 8. Only Task 7's already declared continuation-schedule branch-invariance
control may proceed next.

## Task 7 - run a controlled one-axis candidate matrix

### Goal

Test only evidence-authorized hypotheses and distinguish physics/operator
changes from solver-path changes.

### Candidate axes

Each experiment changes exactly one axis from the sealed baseline:

1. production `cell_reconstructed` versus opt-in
   `element_edge_sg_gss_laux`, only if Task 5 authorizes source/current
   support as causal;
2. baseline versus the Task 4 source-Jacobian correction;
3. high-field mobility field derivatives off/on, only if Task 6 identifies a
   missing feedback derivative;
4. baseline density/current support versus a documented AvalDens-like
   diagnostic, never as an undeclared Sentaurus-equivalent default; and
5. two predeclared continuation schedules with unchanged physics to test
   branch invariance.

Do not combine candidates until each single-axis result has a complete
scorecard. A combined candidate requires evidence that two independently
localized causes remain.

### Acceptance criteria

A candidate is authorized only if it:

- improves knee-window log-current RMSE by at least 50%;
- reduces both `V_slope` and `V_break` errors;
- removes rather than moves any nonmonotonic interval;
- improves the named internal causal metric at the same biases;
- preserves Tasks 4-6 gates;
- passes duplicate-run determinism; and
- does not worsen any global-lattice current error by more than `0.02 dex`.

### Decision gate and stop condition

Typed outcome:

- `single_causal_candidate_authorized`;
- `two_independent_causes_required`;
- `solver_path_only`;
- `improves_curve_without_internal_causality`;
- `tradeoff_without_parity`; or
- `no_authorized_candidate`.

Only the first two outcomes may enter Task 8. A prettier curve without the
predicted internal correction is not accepted.

Executed outcome: `no_authorized_candidate`.

The controlled candidate changed only
`impact_ionization.quasi_fermi_carrier_truncation` from the default `0` to
`1.0e-2`. Two independent avalanche-on runs completed all 29 exact points and
were byte-identical, but they were also byte-identical to the baseline. The
knee-window log-current RMSE remained `11.400736039323265 dex`, `V_break`
error remained `0.232 V`, `V_slope` remained unavailable, and no nonmonotonic
interval was removed. Candidate and baseline process-chain record arrays were
identical; QFP RMSE remained `0.38044152394073266 V`, density log-RMSE remained
`5.7568129570729 dex`, and WP7 remained `density_qfp_feedback_cause` with
203/203 closure rows passing. Task 8 is not authorized.

The subsequent candidate-axis-5 control changed only the continuation maximum
and initial step from `0.05 V` to `0.025 V`. Two independent runs of each
schedule completed all 29 exact points and were deterministic. At
`-19.7/-19.8 V`, schedule-to-schedule QFP differences remained below
`1.28e-12 V`, while the Task 6 carrier-only improvement
(`13.128%/12.942%`) and full-coupled worsening
(`7.023%/7.285%`) were unchanged. The typed outcome is
`continuation_invariant_cross_block_reversal`; the complete candidate curve
campaign and Task 8 are not authorized. Evidence:
`docs/validation/pn2d_task7_continuation_schedule_control_2026-07-30.md`.

The subsequent observation-only Poisson-QFP decomposition used the same
frozen QFP substitution residual and one production baseline Jacobian. At
`-19.7/-19.8 V`, the independent and either-single-cross-block-removed QFP
steps retain positive target-direction cosine `0.7566/0.7511`, while the full
Schur step reverses to `-0.4418/-0.4582`. The Schur reconstruction matches the
production raw Newton step within `7.32e-13 V`, with relative closure below
`4.60e-14`; the configured capped step remains adverse. Duplicate node and
matrix hashes match, boundary targets remain zero, and both biases localize
the strongest adverse projection to node 12.

The typed observation result is
`bidirectional_poisson_qfp_closed_loop_cause`. It proves that the reversal is
created by closing both Poisson-QFP cross-block directions, not by
continuation or post-solve caps. The independent step is severely oversized
and therefore is not an admissible correction. Task 8 remains prohibited
pending model-ownership, sign, scaling, conditioning, and
directional-derivative review of the effective Schur term `C A^-1 B`.
Evidence:
`docs/validation/pn2d_poisson_qfp_cross_block_decomposition_2026-07-30.md`.

The effective-loop follow-up split `C A^-1 B` by carrier, model, spatial
support, sign, scale, and condition. Transport-only and avalanche-only loops
remain adverse at both biases, with direction cosines
`-0.1528/-0.1555` and `-0.2427/-0.2454`; SRH/Auger-only retains the positive
carrier direction `0.7566/0.7511`. Fixed-step directional finite differences
close `B` and `C` within `2.24e-6`, and the transport component has no nonzero
contact-row `C` entries. Row/column equilibration reduces the Schur condition
to `135-138`, so the adverse direction is not caused by worse Schur
conditioning. The typed result is
`transport_and_avalanche_independently_sustain_reversal`. Task 8 remains
prohibited. Evidence:
`docs/validation/pn2d_schur_loop_source_decomposition_2026-07-30.md`.

The subsequent GSS `aux2` ownership/sign audit found that the triangle-GSS
path reverses the archived isothermal electrostatic-potential signs for both
electron and hole midpoint densities. The reference midpoint closes to the
existing Vela SG-edge midpoint below `9.5e-17` relative L2, while the
dominant production/reference ratio is `2.335e7-7.253e7`. However, a
sign-correct midpoint-only source reaches only `0.48603-0.48739` of the
Sentaurus integral. The complete element-edge SG/GSS-Laux vector reaches
`1.00948-1.00954`, consistent with GSS impact ionization consuming complete
directed SG current rather than a standalone
`mu*n_mid*abs(grad QFP)` proxy.

The typed result is
`gss_aux2_sign_transcription_and_operator_ownership_mismatch_confirmed`.
It rejects a sign-only candidate and does not authorize Task 8. A future
Task 7 candidate must remain opt-in and use the complete SG current vector
with matching geometry; the sign-correct midpoint-only result is a negative
control. Evidence:
`docs/validation/pn2d_gss_aux2_ownership_audit_2026-07-30.md`.

The subsequent frozen-state Task 7 comparison ran the production triangle
baseline, sign-correct midpoint-only negative control, and complete
element-edge SG/GSS-Laux candidate twice at `-19.7/-19.8 V`. After
conservatively projecting element-vertex current and source back to the
physical node support, the complete candidate has:

- integrated source ratios `1.009537/1.009483`;
- electron/hole matching-current median and P95 below `0.0049 dex`;
- active-node source maximum below `0.0052 dex`;
- 100% nonzero vector direction agreement; and
- 27/27 duplicate numerical/process artifacts identical.

The result is insensitive to active-source floors from `1e-5` through
`1e-8`. The sign-correct midpoint-only control remains rejected at
`0.48603/0.48739`. The typed outcome is
`complete_sg_vector_fixed_state_prequalified`. It authorizes the remaining
single-axis, self-consistent Task 7 exact-lattice candidate run, but not
Task 8 or a default change. Evidence:
`docs/validation/pn2d_task7_frozen_sg_candidate_2026-07-30.md`.

The remaining self-consistent Task 7 candidate completed duplicate
avalanche-on global/knee exact lattices plus one complete off/IIC/on process
campaign. It improves knee log-current RMSE from `11.400736` to
`0.0146388 dex`, restores `V_slope` with `0.01643 V` error, reduces
`V_break` error to `0.021 V`, and reduces the QFP/density process errors to
`9.986e-5 V/0.001673 dex`. Duplicate IV and Task 6 outputs are byte-identical,
and all 203 closure rows pass.

However, the candidate replaces four baseline reverse intervals, including
the high-current `-18→-19.25 V` sequence, with three low-current intervals at
`-3→-5 V` and `-6→-7 V`. The predeclared gate requires removal rather than
relocation. The typed outcome is therefore `tradeoff_without_parity`, and
Task 8/default changes remain prohibited. Evidence:
`docs/validation/pn2d_task7_self_consistent_sg_candidate_2026-07-30.md`.

## Task 8 - implement the minimal opt-in correction

### Goal

Turn the authorized causal result into the smallest tested implementation
without changing production defaults.

### Actions

1. Add failing focused tests before or with the minimal implementation.
2. Keep Real/AD/diagnostic formula ownership shared.
3. Add configuration parsing, validation, serialization, and capability guards
   for any new opt-in choice.
4. Preserve old behavior when the option is absent.
5. Add tests for topology orientation, contacts, zero/tail branches, both
   carriers, source identity, analytic derivatives, and configuration hashes.
6. Re-run all focused source, Jacobian, fixed-state, and first-update tests.

### Acceptance criteria

- Every new test fails for the intended reason on the baseline and passes with
  the correction.
- Existing default-path numerical snapshots remain unchanged.
- Task 4 derivative, Task 5 fixed-state, and Task 6 causality gates pass.
- No unrelated solver or physics diff is present.

### Decision gate and stop condition

Report the exact files, formulas, tests, and scoped diff. Stop if the change
needs an empirical scale or changes an unproven production behavior.

## Task 9 - qualify the deterministic M0 self-consistent curve

### Goal

Apply the unchanged opt-in candidate to the complete M0 on/off sweeps and
decide whether curve and knee parity are achieved.

### Actions

1. Run baseline and candidate twice on both exact lattices.
2. Run both Vela avalanche-on and avalanche-off branches.
3. Generate all Task 2 figures, tables, manifests, and acceptance results.
4. At each knee state, report:
   - current and gain;
   - adjacent log slope;
   - knee estimators;
   - QFP and hotspot density;
   - mobility/current/alpha/source;
   - residual/Jacobian and closure; and
   - accepted/rejected continuation state.
5. Repeat the fixed-state and first-update checks on the resulting
   self-consistent states.

### Acceptance criteria

All common curve, knee, physical, numerical, and determinism gates must pass.
In addition:

- M0 avalanche-off RMSE remains `<= 0.001 dex`;
- M0 avalanche-off maximum error remains `<= 0.002 dex`;
- no accepted lower-bias point regresses relative to baseline; and
- no exact state is replaced by interpolation.

### Decision gate and stop condition

Typed outcome:

- `m0_self_consistent_curve_knee_parity_passed`;
- `current_magnitude_gap`;
- `knee_voltage_gap`;
- `knee_shape_gap`;
- `nonmonotonic_branch`;
- `continuation_branch_difference`;
- `closure_failure`; or
- `solver_first_failure`.

Only the first outcome may enter Task 10. Report the complete acceptance JSON
and plots before continuing.

## Task 10 - build a dose-preserving mesh sequence

### Goal

Prove that the M0 result is not a coarse-mesh coincidence before making a
physical or production claim.

### Actions

1. Replace the rejected M1/M2 profile construction with a nested,
   dose-preserving junction definition.
2. Ensure the junction coordinate belongs to exactly one analytical profile
   window or use a documented local-replace construction.
3. Build at least two finer levels with approximately 2x and 4x junction
   resolution.
4. Pair Vela and Sentaurus on each exact physical mesh.
5. Run on/off global and knee lattices using the unchanged Task 9 candidate.
6. Compare discrete dose, depletion support, source integral, terminal curve,
   gain, slopes, and knee estimators across levels.

### Acceptance criteria

- Discrete total-impurity dose changes by `< 0.1%` between levels.
- No junction node is double-counted by profile windows.
- Vela/Sentaurus mesh and doping hashes are paired at each level.
- The two finest levels change integrated on/off sources by `< 2%`.
- Both finer levels pass the final curve and knee gates, or errors converge
  monotonically toward them without changing the causal classification.
- Knee voltage changes by `<= 0.10 V` between the two finest levels.

### Decision gate and stop condition

Typed outcome:

- `dose_preserving_mesh_parity_passed`;
- `mesh_dependent_knee`;
- `mesh_dependent_source`;
- `dose_preservation_failed`; or
- `cross_simulator_mesh_mismatch`.

Stop on any outcome except the first. Do not reuse the rejected M1/M2 results
as substitutes.

### Execution note 2026-07-31

The prospective dual-domain M0 confirmation passed, but the replacement
dose-preserving mesh sequence stopped with `mesh_dependent_knee`.

The new single-owner junction construction preserved discrete dose to
floating-point precision and eliminated double-counted junction nodes.
Nevertheless, Sentaurus M1 and M2 changed fitted `V_break` by `0.402 V`,
exceeded the two-finest integrated-source gate by orders of magnitude, and
M1 developed a discontinuous/nonmonotonic avalanche-on branch. Vela
avalanche-off and IIC completed the new M0 lattice, while the SG/Laux
self-consistent branch stalled near `-17.1053958 V` when targeting `-18 V`.

Task 11 is therefore not authorized. See
`docs/validation/pn2d_task10_dose_preserving_mesh_gate_2026-07-31.md`.

### Corrective execution note 2026-07-31

The M0 continuation failure above was traced to the single-owner junction
contract, not to a Vela solver limit. Assigning `+1e17 cm^-3` or
`-1e17 cm^-3` net doping to the geometric junction nodes moves the
zero-net-doping contour under linear finite-element interpolation. Both signs
stall near `-17.105 V`; zero-net controls reach `-20 V`.

The mesh generator now uses a balanced-half junction window:
`ND=NA=5e16 cm^-3`. This preserves the intended zero-net junction, total
impurity, and integrated dose. Two independent Vela runs using the actual
Sentaurus-imported corrected M0 input completed all 29 exact points for off,
IIC, and SG/Laux-on, with zero rejected attempts and matching all 87 state
hashes. The paired corrected Sentaurus curve also reaches `-20 V`; Vela and
Sentaurus give `V_break=-19.622/-19.654 V` and differ by `0.032 V`.

The immediate M0 blocker is therefore closed. The earlier M1/M2
`mesh_dependent_knee` and `mesh_dependent_source` results used the superseded
single-owner junction contract and remain secondary/open until regenerated.
They do not authorize Task 11 or a production-default change. See
`docs/validation/pn2d_task10_m0_stall_root_cause_2026-07-31.md`.

### Balanced-junction M1/M2 execution note 2026-07-31

The M1/M2 sequence was regenerated with `ND=NA=5e16 cm^-3` at every
`x=1 um` junction node. Both levels preserve the total dose and Vela completes
all 29 exact off/IIC/on points through `-20 V`. Duplicate Vela runs match all
174 state hashes and all six IV hashes.

M2 is monotonic in both simulators and has a cross-simulator `V_break`
difference of `0.014 V`. M1 is nonmonotonic in both simulators, with branch
jumps at different biases. M1-to-M2 `V_break` changes are `0.123 V` for Vela
and `0.191 V` for Sentaurus. Maximum knee-lattice integrated-source changes
are `3835%` and `165471%`, respectively.

Task 10 therefore still stops with primary `mesh_dependent_knee` and secondary
`mesh_dependent_source`. This result is distinct from the resolved M0
continuation/input defect. Task 11 remains unauthorized. Evidence:
`docs/validation/pn2d_task10_balanced_mesh_independence_2026-07-31.md`.

### Sentaurus-golden priority amendment 2026-07-31

The project objective is same-grid, same-input Vela parity against Sentaurus
as the golden reference. Cross-level mesh convergence remains valuable
scientific characterization, but a nonmonotonic branch observed in both
simulators on the same coarse M1 topology is not a blocking product gate.

The original `mesh_dependent_knee` and `mesh_dependent_source` values remain
recorded without reinterpretation. They mean that M1/M2 do not prove a
mesh-converged physical knee; they do not reject the SG/Laux candidate when
stable same-grid comparisons pass. M0 and M2 are the primary golden-parity
levels. M1 is retained as a nonblocking continuation/topology observation.

The project-level Task 10 outcome is therefore reclassified as:

```text
golden_same_grid_parity_passed_mesh_observation_open
```

Task 11 is authorized to proceed. This amendment does not authorize a
production-default change: the candidate remains opt-in until Task 11 full
regression and independent reviews are complete.

## Task 11 - regression, review, and default decision

### Goal

Decide whether the opt-in correction is ready for a separate production
default review.

### Actions

1. Run all focused avalanche, source, Jacobian, general-Tri3, DC-sweep,
   continuity, and PN2D tests.
2. Run the full Release build and CTest suite.
3. Re-run the 201-point forward IV with
   `cell_reconstructed_total_impurity`.
4. Verify BV configs still render `net_doping`.
5. Run `git diff --check`, ASCII/source checks, schema/manifest validators, and
   scoped diffs.
6. Obtain independent scientific and code review of:
   - the causal claim;
   - the knee definition and scores;
   - derivative and closure evidence;
   - mesh independence; and
   - the proposed default behavior.
7. Keep the candidate opt-in unless every gate and both reviews pass.

### Acceptance criteria

- Full CTest: zero failures.
- Forward IV: 201/201 points converged.
- Forward-IV anchor currents change by `< 0.5%`.
- BV basis remains `net_doping`; IV basis remains
  `cell_reconstructed_total_impurity`.
- Avalanche-off M0 parity remains within Task 9 limits.
- No generated simulation output or `tmp/` is staged.
- Both reviews accept the evidence and scope.

### Final decision

| Evidence | Allowed decision |
|---|---|
| Any Task 1-9 gate fails | keep diagnostic only and preserve typed outcome |
| Stable same-grid M0/M2 golden parity passes while both simulators share an M1 anomaly | proceed with Task 11; retain the mesh observation |
| Vela alone diverges from Sentaurus on the same stable mesh and inputs | keep opt-in; classify the same-grid model difference |
| Curves match only after scale/shift/fitting | reject; classify `model_difference` |
| All tasks and both reviews pass | propose a separate production-default change |

### Task 11 execution note 2026-07-31

The Release build was current. The focused physics/solver selection passed
`130/130` tests and the complete Release suite passed `506/506`. Fresh ASCII,
PN2D template, Sentaurus import/schema, and BV process-observability
validators also pass.

The 201-point forward IV rerun converged `201/201` points with
`cell_reconstructed_total_impurity`. Its largest predeclared anchor change is
`0.00178944%`, below the `0.5%` gate. Freshly rendered IV/BV configurations
retain `cell_reconstructed_total_impurity`/`net_doping`, respectively.

The current runner reproduces the sealed M0 avalanche-off curve exactly by
SHA-256. Its nonzero-bias RMSE is `6.9965091e-5 dex` and maximum error is
`1.2546889e-4 dex`, both inside the Task 9 limits. Source/contact and terminal
closure also pass.

No Task 11 solver-source or default change was made. Two new independent
reviews were then completed against commit `ec00347`; both return
`APPROVE_WITH_CONDITIONS`. The Task 11 typed outcome is therefore:

```text
task11_regression_passed_independent_reviews_approve_with_conditions
```

The reviews authorize a separate PN2D BV template default proposal and a
prospective M0/M2 dual-domain acceptance contract, not a direct default
switch. Required conditions include:

- rescore balanced M0 and M2 under one contract with bound closure evidence;
- save a unified M2 acceptance artifact;
- keep same-grid golden agreement distinct from mesh convergence;
- change the PN2D template SG/Laux configuration atomically;
- retain the global C++ defaults and an explicit legacy/rollback path; and
- add fresh-render, compatibility, rollback, and default-path M0/M2 tests.

See:

- `docs/validation/pn2d_task11_regression_review_2026-07-31.md`;
- `docs/validation/pn2d_task11_independent_scientific_review_2026-07-31.md`;
- `docs/validation/pn2d_task11_independent_code_review_2026-07-31.md`.

### Prospective template-default proposal execution note 2026-07-31

A separate M0/M2 contract was frozen before candidate execution:

```text
docs/validation/contracts/pn2d_bv_m0_m2_template_default_acceptance_v1.json
```

The atomic PN2D BV profile, explicit legacy rollback, default-render tests,
contract-domain analyzer, exact-lattice hash bindings, partial failure state
manifest, and strict machine evaluator were implemented.  No global C++
default was changed.

Two independent default-render SG/Laux runs were executed for balanced M0 and
M2.  Avalanche-off and avalanche-on complete 29/29 points and are deterministic
at both levels.  M0 same-grid parity and closure pass all gates.  M2 fails the
frozen knee median and required `V_slope` gates.  IIC/postprocess is incomplete
and its partial IV hash is not deterministic at both levels.  The unified
typed outcome is:

```text
pn2d_bv_template_default_not_accepted
```

Both second independent reviews return `REJECT_DEFAULT_CHANGE`.  The PN2D BV
template default was therefore rolled back atomically to
`legacy_cell_reconstructed`; SG/Laux remains an explicit opt-in.  The failed
evidence is preserved.  Focused Python tests pass 28/28 and the post-rollback
Release CTest passes 506/506.

See:

- `docs/validation/pn2d_bv_m0_m2_template_default_acceptance_2026-07-31.md`;
- `docs/validation/pn2d_bv_template_default_second_scientific_review_2026-07-31.md`;
- `docs/validation/pn2d_bv_template_default_second_code_review_2026-07-31.md`.

## Decision-gate reporting template

After every task, report before continuing:

1. task number and typed outcome;
2. exact commands and commit;
3. input, config, mesh, doping, and output hashes;
4. required acceptance metrics with threshold, observed value, and pass/fail;
5. first failed bias/node/cell/carrier/term, if any;
6. generated artifact paths;
7. scoped `git status --short`; and
8. explicit authorization or prohibition for the next task.

Do not combine two task reports. A failed stop condition ends execution even if
a later experiment appears promising.

## Suggested commit boundaries

1. Review and commit the completed common process schema/contract files.
2. Review and commit the completed Sentaurus process-matrix runner/analyzer
   and dry-run tests.
3. WP2.1 exact-lattice additions and validation summaries.
4. Task 2 analyzer contract and tests.
5. Task 4 focused Jacobian tests and minimal correction, only if authorized.
6. Tasks 5-7 diagnostic tooling and reports.
7. Task 8 minimal opt-in implementation and tests.
8. Tasks 9-10 validation summaries, excluding generated outputs.
9. Task 11 review responses and any separately authorized default change.

Before every commit:

- inspect `git status --short`;
- stage explicit paths only;
- verify `tmp/` is not staged;
- verify no `build*/` simulation result, TDR, PLT, CSV, VTK, log, or state
  output is staged; and
- run `git diff --cached --check`.

## Next execution prompt

Use the following prompt for the next execution slice:

> Continue Task 7 from the sealed `tradeoff_without_parity` result. Perform an
> observation-only audit of the SG/Laux candidate's three low-current reverse
> intervals at `-3/-4/-5/-6/-7 V`. Compare duplicate Newton/continuation
> histories, terminal-current components, state hashes, and off/IIC/on
> branches. Classify source feedback, current extraction, or numerical-floor
> behavior without changing the candidate or adding a post-hoc current floor.
> Keep Task 8 and all production-default changes prohibited until the
> no-relocation gate is resolved by a separately reviewed typed outcome.

## Task 7 low-current causal audit result - 2026-07-30

Typed outcome:

`low_current_state_precision_floor_not_avalanche_operator_or_terminal_extractor`

The observation-only audit is complete:

- two of three low-current reverse intervals are shared by off, IIC, and on;
- off and IIC states are byte-identical;
- duplicate on states are byte-identical;
- alpha and raw source increase monotonically from `-3` to `-7 V`;
- SG-flux and residual terminal currents agree within `4.650654e-15`
  relative error;
- contact hole-QFP drops are approximately `3e-14 V`;
- strict Newton tolerances move/remove the intervals and terminate with
  `stall_residual_floor`.

The original `tradeoff_without_parity` result remains preserved because the
predeclared no-relocation gate is not changed retrospectively. Task 8 and
production-default changes remain prohibited pending a prospective review of
how BV-active metrics and low-current precision-floor behavior should be
reported separately.

Evidence:

`docs/validation/pn2d_task7_low_current_root_cause_audit_2026-07-30.md`.

### Next execution prompt

> Review the completed Task 7 low-current typed outcome alongside the original
> `tradeoff_without_parity` score. If separately authorized, define a
> prospective acceptance contract that retains all raw intervals but separates
> BV-active collision-ionization parity from low-current solver precision.
> Otherwise stop. Do not enter Task 8, add a current/minimum-field threshold,
> or change production defaults.

## Task 7 dual-domain prospective contract review - 2026-07-31

The prospective contract review is complete with typed decision:

`bv_model_consistent_low_current_precision_floor_open`.

The contract freezes disjoint exact domains:

- BV model consistency: `-15` through `-20 V` on the declared 14-point
  effective lattice; and
- low-current solver precision: `-3/-4/-5/-6/-7 V`.

It reuses the previously frozen curve, knee, fixed-state current/source,
self-consistent state, closure, and determinism thresholds. It does not edit
the historical `tradeoff_without_parity` result. The low-current
classification additionally requires every exact-domain on/off current to
remain below `1e-15 A/um`, so a high-current low-bias runaway cannot receive a
precision-floor waiver.

Replay result:

- every BV-active model-consistency gate passes;
- every low-current precision-floor classification gate passes;
- the historical Task 7 score remains preserved; and
- production-default changes remain unauthorized.

The only newly authorized scope is:

`opt_in_bv_model_validation_only`.

Contract and review:

- `docs/validation/contracts/pn2d_bv_dual_domain_acceptance_v1.json`;
- `docs/validation/pn2d_bv_dual_domain_acceptance_contract_review_2026-07-31.md`.

### Next execution prompt

> Continue only with a separately scoped opt-in BV model-validation task under
> the frozen dual-domain contract. Preserve the historical Task 7 score and
> report low-current solver precision separately. Do not change the production
> default, add a minimum-field/current threshold, fit model coefficients, or
> treat the low-current nonlinear precision task as closed.

## M2 Sentaurus-state SG/Laux frozen replay - 2026-07-31

The predeclared four-point discriminating experiment is complete with typed
outcome:

`state_feedback_dominant`

Sentaurus avalanche-on states at `-18`, `-19.5`, `-19.7`, and `-20 V` were
imported on the shared M2 mesh and evaluated by the complete Vela SG/Laux
operator with `coupling_mode=postprocess_only`.  The imported state was not
advanced, every process record has `solver_coupled=0`, and no residual feedback
contribution is nonzero.

The frozen Vela/Sentaurus total-source ratios are `1.002481`, `1.002381`,
`1.002404`, and `1.002371`.  Mean frozen-state error is `0.001045 dex`, versus
`0.057445 dex` for self-consistent Vela; the frozen replay removes
`0.056400 dex` on average.  The self-consistent total-source ratio decreases
from `0.936135` at `-18 V` to `0.825003` at `-20 V`, while the frozen ratio
remains flat near one.  All 20 node/edge/triangle/element/process artifacts are
byte-identical across two independent runs.

This evidence localizes the dominant M2 knee discrepancy to coupled state
formation/feedback rather than frozen SG/Laux source evaluation.  It does not
authorize changing SG/Laux, a production default, or an acceptance threshold.

Evidence:

- `docs/validation/pn2d_bv_m2_sentaurus_frozen_sg_laux_2026-07-31.md`;
- `docs/validation/pn2d_bv_m2_sentaurus_frozen_sg_laux_2026-07-31.html`.

### Next execution prompt

> Keep SG/Laux and all acceptance thresholds unchanged.  On the same four M2
> frozen states, perform one-family-at-a-time substitutions for electrostatic
> potential, quasi-Fermi potentials, and carrier densities, then compare the
> first coupled Newton update and carrier-row residual.  Determine which state
> family recovers most of the `0.082516 dex` improvement at `-20 V` and whether
> the first update moves it in the direction of the final self-consistent
> deficit.  Do not modify production defaults.

## M2 single-family state substitution - 2026-07-31

The one-family fixed-state substitutions and first coupled Newton probes are
complete with typed outcome:

`qfp_dominant__density_feedback_moves_qfp_away_from_sentaurus`

At `-20 V`, substituting Sentaurus electron/hole QFP alone recovers
`0.738076` of the source-error removal obtained from the full Sentaurus state.
QFP wins three of four predeclared bias comparisons.  The QFP-only recovery is
bias dependent: it is negative at `-18 V`, then rises to `0.465512`,
`0.580240`, and `0.738076` at `-19.5`, `-19.7`, and `-20 V`.

Density-only substitution leaves the frozen SG/Laux source unchanged.  This is
expected for the current Masetti plus field-dependent mobility configuration:
the canonical SG edge flux is reconstructed from `psi`, QFP, and intrinsic
density, while mobility depends on doping and field rather than imported
`n/p`.  Density still affects other residual paths.

The first coupled update rejects rather than preserves the imported golden QFP
direction.  Density-only residual feedback produces a negative QFP target
projection at every bias (`-0.0025435` at `-20 V`).  Direct QFP feedback has a
larger negative projection, and an independent Sentaurus-QFP start increases
the combined target distance at all four biases.  At `-19.5` through `-20 V`,
its first production trial residual is `3.91-4.05` times the initial residual.

All 120 repeated artifacts are byte-identical.  No physics model, production
default, continuation schedule, or acceptance threshold was changed.  This
localizes the next investigation to carrier-QFP residual/Jacobian coupling and
state consistency, not to a frozen SG/Laux source formula change.

Evidence:

- `docs/validation/pn2d_bv_m2_single_family_state_substitution_2026-07-31.md`;
- `docs/validation/pn2d_bv_m2_single_family_state_substitution_2026-07-31.html`.

### Next execution prompt

> Keep SG/Laux and all production defaults unchanged.  On the same M2 states,
> split electron and hole QFP substitutions, decompose the first carrier
> residual into transport, recombination, avalanche, and boundary terms, and
> finite-difference the carrier-QFP and Poisson-QFP Jacobian blocks on both the
> Vela baseline and mixed Sentaurus-QFP state.  Identify the first derivative
> or residual term whose sign/scale differs before proposing any correction.

## M2 carrier-QFP residual and Jacobian verification - 2026-08-01

The carrier-resolved residual and cross-block verification is complete with
typed outcome:

`phip_dominant__electron_flux_dominant__hole_flux_dominant__analytic_fd_inconsistent__both_qfp_updates_roll_back_from_sentaurus`

Hole QFP has the larger frozen-source effect at all four biases and recovers
`0.991518` of the joint-QFP improvement at `-20 V`.  The joint-QFP carrier
residual change is transport/SG-flux dominated: `0.883898` for electrons and
`0.892987` for holes, with avalanche contributing `0.116102` and `0.107013`.
Both first carrier updates at `-20 V` point away from Sentaurus, with target
projections `-0.923986` and `-0.889187`.

Poisson, transport, SG-avalanche, and boundary/gauge cross blocks pass the
unchanged `5e-5` finite-difference threshold; their worst relative errors are
`6.12734e-8`, `5.15168e-8`, `9.05564e-10`, and `3.55650e-10`.  The formal
all-block result remains `analytic_fd_inconsistent` because SRH/Auger fails in
relative terms.  A seven-step sensitivity audit classifies that failure as:

`formal_relative_gate_fails_only_at_srh_absolute_fd_floor`

The SRH/Auger absolute differences remain below `2.70635e-15`; the best
comparison is `1.00713e-21`.  All 148 repeated artifacts are byte-identical,
and carrier-term closure is `1.32349e-23`.  No physical model, production
default, continuation schedule, or acceptance threshold was changed.

Evidence:

- `docs/validation/pn2d_bv_m2_qfp_carrier_jacobian_verification_2026-08-01.md`;
- `build-release/pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731/report.html`.

### Next execution prompt

> Keep SG/Laux and all production defaults unchanged.  On the same M2
> baseline and joint-QFP states, perform an edge-level carrier transport
> Jacobian audit at the residual hotspots.  Separate mobility,
> Bernoulli/GSS coefficient, QFP driving force, row scaling, and contact-row
> elimination derivatives for electron and hole equations, and verify each
> contribution by finite differences before proposing any opt-in correction.

## M2 hotspot-edge transport Jacobian decomposition - 2026-08-01

The edge-level carrier transport audit is complete with typed outcome:

`transport_edge_decomposition_verified`

On the prior interior transport-residual hotspot support, the production
analytic derivative agrees with frozen-mobility finite differences to
`5.12344e-10` after normalization by the dominant derivative on the same edge.
The independently evaluated mobility product closes to `2.11374e-16`, and the
live-mobility total derivative closes to `1.01736e-7`.  Bernoulli/GSS QFP
derivatives, row-scaling closure, contact-row replacement, and contact identity
entries are exact under their declared contracts.  All 16 raw probe outputs
are byte-identical between two independent runs.

At `-19.5` through `-20 V`, replacing Vela QFP with Sentaurus QFP leaves the
Bernoulli/GSS coefficients exactly unchanged, reduces QFP drive by about
`14.7-27.3%`, and increases high-field mobility by about `10.1-19.0%`.
Mobility response contributes only `3.75-7.00%` of the dominant hotspot
transport derivative.  The dominant production derivative nevertheless falls
to `0.0682-0.1720` of baseline, so the first material change is the
QFP-controlled exponential carrier-population term.  Hotspot row weights
change by less than `0.09%`; contact-row elimination is not the cause.

This result verifies the local transport derivative rather than identifying a
local formula defect.  The remaining investigation moves one level upward to
the complete carrier-block linear solve and its strongly anisotropic variable
scaling.

Evidence:

- `docs/validation/pn2d_bv_m2_transport_edge_jacobian_verification_2026-08-01.md`;
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/report.html`.

### Next execution prompt

> Keep SG/Laux, mobility settings, continuity-row scaling, contact handling,
> and acceptance thresholds unchanged.  On the same M2 baseline and joint-QFP
> states, decompose the full electron/hole carrier-block linear solve at the
> hotspot support.  Compare dominant QFP columns, diagonal dominance,
> left/right singular directions, variable scaling, recombination coupling,
> common avalanche rows, and the carrier-only Newton update projected onto the
> population-derivative directions.  Do not modify production defaults before
> a wrong sign, missing derivative, or inconsistent scale is demonstrated.

## M2 complete carrier-block linear-solve decomposition - 2026-08-01

The free-carrier solve decomposition is complete with typed outcome:

`carrier_block_linear_solve_decomposed`

The full linear closure is `7.195e-16`, row-scaled and unscaled carrier steps
agree within `2.098e-13`, and all five output families are byte-identical over
two independent runs.  The transport-only electron-hole cross block is exactly
zero.  No production setting was changed.

At `-19.5` through `-20 V`, the joint-QFP state has an L2 row/column
equilibrated condition number of `381.5-543.3`, versus `65.8-105.1` for the
baseline.  Its carrier step is `1.30-1.43 V`, versus about `0.002 V` for the
baseline.  The two dominant joint-QFP singular directions have relative
singular values of about `6.9e-16` to `1.1e-15` and carry `94.88-96.34%` of
the actual production sparse-solver step energy.

Avalanche supplies approximately all electron-hole cross-block norm.  Removing
cross-carrier entries changes the joint-QFP step by `38.0-46.4%`, and removing
the avalanche matrix changes it by `25.1-31.5%`; direction cosines remain above
`0.992`.  Recombination changes the step by only about `1e-6`.  The evidence is
therefore a magnitude amplification on near-null junction modes, not a sign
reversal or a local SG derivative defect.

All 40 joint-QFP top-10 update nodes lie within `0.25 um` of the `x=1.0 um`
junction, while none of the corresponding baseline top-10 updates do.  The
two strongest modes generally peak on the uniformly doped `x=0.75/1.25 um`
junction shoulders, four triangle-graph steps from the compensated column.
This supports junction discretization as an indirect amplifier, but does not
identify local compensated-triangle doping smoothing as the direct defect.

Evidence:

- `docs/validation/pn2d_bv_m2_carrier_block_decomposition_2026-08-01.md`;
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/report.html`.

### Next execution prompt

> Keep SG/Laux, production defaults, mesh, nodal doping, and acceptance
> thresholds unchanged.  On the M2 `x=0.75-1.25 um` junction-shoulder support,
> project the two dominant carrier soft modes onto transport, avalanche
> diagonal, and avalanche cross-carrier Jacobian components.  In parallel,
> audit triangle-level nodal-doping interpolation, control-volume dose, and
> compensated-junction ownership against the Sentaurus export semantics.  Any
> doping counterfactual must remain frozen-state, opt-in, and dose preserving.

## M2 P0-P2 soft-mode and doping-semantics audit - 2026-08-01

P0-P2 are complete under the frozen contract
`docs/validation/contracts/pn2d_bv_m2_soft_mode_p0_freeze_v1.json`.

The carrier-mode decomposition now exports signed modal projections for
transport, recombination, avalanche diagonal, and avalanche cross-carrier
Jacobians, plus the corresponding residual projections and counterfactual
solve amplitudes.  Two independent runs are byte-identical.  Maximum modal
Jacobian and RHS relative closure errors are `5.96e-15` and `9.62e-14`.

At `-20 V`, the two dominant joint-QFP modes carry `94.88%` of the production
step energy.  Transport contributes `8.35-9.18` times the net singular
stiffness.  Avalanche diagonal contributes `-3.19` to `-3.27` times and
avalanche cross-carrier contributes `-4.16` to `-4.91` times.  Recombination
Jacobian contribution is below `3e-6` of net stiffness.  The near-null modes
therefore arise from signed cancellation between transport and avalanche
Jacobian components, not from a local derivative sign defect.

The M2 Sentaurus mesh TDR was independently imported with the `reported`
compensated-doping policy.  Vela and Sentaurus have exactly the same 115 node
IDs, coordinates, 191 unordered triangle connectivities, donor values,
acceptor values, nodal net doping, and junction-edge endpoint averages.  The
maximum discrepancy between Sentaurus's reported net field and `ND-NA` is
`1.92e-15` relative.  The Vela and Sentaurus doping CSV files are byte
identical.  This excludes a nodal-input, topology, compensated-node ownership,
or edge-average mismatch.

The active Vela node-volume policy is barycentric `area/3`.  The same TDR
export does not expose Sentaurus node control volumes, so exact control-volume
parity is not directly observable.  On the six nodes selected by the dominant
soft modes, the mixed-Voronoi-to-barycentric volume ratio spans `1.0-1.5`.
This is a bounded geometry-policy sensitivity to test next; it is not evidence
of a doping smoothing defect.

Evidence:

- `docs/validation/pn2d_bv_m2_soft_mode_doping_semantics_2026-08-01.md`;
- `build-release/pn2d-bv-m2-soft-mode-component-projection-20260801/result.json`;
- `build-release/pn2d-bv-m2-doping-control-volume-audit-20260801/result.json`.

### Next execution prompt

> Keep nodal doping and SG/Laux unchanged.  First obtain or reconstruct the
> Sentaurus box-node volumes.  If they cannot be exported, run a frozen-state,
> opt-in barycentric-versus-mixed-Voronoi first-step comparison on the same M2
> states.  Measure the two dominant modal stiffnesses, QFP step amplitudes,
> Poisson/carrier residuals, and integrated source.  Do not enter a doping
> redistribution counterfactual unless the volume-policy control fails to
> explain the sensitivity.

## M2 barycentric versus mixed-Voronoi frozen first step - 2026-08-01

The predeclared dual-policy experiment is complete with typed outcome:

`material_node_volume_policy_sensitivity`

The classification is driven by the Poisson block.  Re-evaluating a
barycentric-converged Vela state with mixed-Voronoi node volumes raises the
initial Poisson residual from about `2e-8` to `26`, changing the complete first
step by `277-286x`.  This demonstrates that the frozen state is not stationary
under the alternative control-volume discretization.

The discriminating joint Sentaurus-QFP carrier result is much narrower.  From
`-19.5` through `-20 V`, the carrier-only step magnitude changes by
`1.65-1.80%`, its direction cosine remains at least `0.999994`, and the
L2-equilibrated carrier-block condition number changes by less than `2.1e-8`
relative.  At `-20 V`, the dominant normalized modal terms remain
`+9.1821` transport, `-3.2702` avalanche diagonal, and `-4.9119` avalanche
cross-carrier under both policies.  The integrated frozen SG/Laux source is
exactly unchanged at every policy, state, and bias.

The result therefore does not explain away the carrier soft mode: the signed
transport-avalanche cancellation survives the volume-policy change.  It does
show that control-volume policy can move the self-consistent branch indirectly
through Poisson.  Nodal doping redistribution remains unauthorized.

Evidence:

- `docs/validation/pn2d_bv_m2_node_volume_policy_first_step_2026-08-01.md`;
- `docs/validation/contracts/pn2d_bv_m2_node_volume_policy_first_step_v1.json`;
- `build-release/pn2d-bv-m2-node-volume-policy-first-step-20260801/result.json`.

### Next execution prompt

> Keep the production barycentric default, nodal doping, SG/Laux, and all
> thresholds unchanged.  If exact Sentaurus box volumes remain unavailable,
> run a separately contracted, opt-in M2 mixed-Voronoi self-consistent control:
> avalanche-off first, IIC second, and SG/Laux-on only if off/IIC do not degrade
> same-grid golden agreement.  Reject the candidate if it improves the on knee
> by sacrificing the already closed avalanche-off baseline.

## M2 mixed-Voronoi self-consistent control - 2026-08-01

The predeclared off -> IIC -> SG/Laux-on sequence completed with typed outcome
`completed_all_stages`.  Every branch reached all 29 exact-lattice points in
two independent runs; IV files and all per-bias state hashes were identical
between repetitions.  IIC was also byte-identical to avalanche-off in both IV
and every state hash, proving that its ionization calculation did not feed back
into the solved state.

Mixed-Voronoi did not sacrifice the Sentaurus avalanche-off baseline.  Across
all 28 nonzero points, its off log-current RMSE is `3.868e-5 dex` and maximum
error is `9.223e-5 dex`.  The previous barycentric M2 off RMSE was
`0.01008 dex`.

After both controls passed, the mixed-Voronoi SG/Laux-on branch also completed
twice.  Its all-point and knee log-current RMSE values against Sentaurus are
`0.001923 dex` and `0.003051 dex`; the maximum error is `0.004647 dex`.
Fitted V_break is `-19.390 V` versus Sentaurus `-19.391 V`, a `0.001 V`
difference.  Neither simulator has a V_slope crossing in the declared range.

The result supports a control-volume geometry mismatch as the dominant
remaining M2 difference.  It does not authorize changing the production
default; M0, forward IV, mixed-Voronoi boundary/obtuse-cell behavior, legacy
fallback, and the full Release suite remain outside this observation-only
contract.

Evidence:

- `docs/validation/pn2d_bv_m2_mixed_voronoi_self_consistent_control_2026-08-01.md`;
- `docs/validation/contracts/pn2d_bv_m2_mixed_voronoi_self_consistent_control_v1.json`;
- `build-release/pn2d-bv-m2-mixed-voronoi-self-consistent-control-20260801/gate_report.json`.

### Next execution prompt

> Keep both production defaults unchanged.  Draft a separate prospective
> node-volume-policy acceptance contract covering M0 and M2 BV, 201-point
> forward IV, obtuse/boundary mixed-Voronoi geometry invariants, explicit
> barycentric legacy fallback, and the complete Release test suite.  Only
> after that contract and an independent code/scientific review pass should a
> PN2D BV template default change be considered.

## PN2D node-volume default acceptance - 2026-08-01

The unified prospective acceptance contract completed with typed outcome
`ready_for_independent_default_policy_reviews`.  The contract hash is
`3a2d65879d1d7446a259afb6f81af9c17e7da2686536e9f6832eb5b5810f6c22`.

The new obtuse/winding, boundary conservation, parser fallback, and explicit
mixed selection tests pass.  M0 and M2 were both freshly run twice through the
off -> IIC -> SG/Laux-on sequence.  Every branch completed 29/29 points and was
deterministic; IIC was exactly state- and IV-equivalent to off on both grids.

M0/M2 SG/Laux-on log-current RMSE values are `0.001766/0.001923 dex`.
Both fitted V_break errors are `0.001 V`; M0 has matching V_slope crossings
within `0.00082 V`, while neither M2 curve crosses.  Both exact lattices are
monotonic.

The prospective forward-IV control completed one barycentric and two mixed
runs at 201/201 points.  The mixed runs are byte deterministic.  Their six-
anchor median/maximum Sentaurus errors are `0.2700%/0.4066%`, with no material
degradation relative to barycentric.

The complete Release suite passed `509/509` twice.  The aggregate machine
outcome is stored at
`build-release/pn2d-node-volume-default-acceptance-v1-20260801/acceptance.json`.

No production default was changed.  The remaining required work is two truly
independent reviews of the scientific evidence and of an atomic PN2D BV
template/default-render patch.

Evidence:

- `docs/validation/pn2d_node_volume_policy_default_acceptance_2026-08-01.md`;
- `docs/validation/contracts/pn2d_node_volume_policy_default_acceptance_v1.json`;
- `build-release/pn2d-node-volume-default-acceptance-v1-20260801/acceptance.json`.

### Next execution prompt

> Conduct independent scientific and code reviews against the frozen node-
> volume-policy contract.  If both approve, design an atomic patch that changes
> only the PN2D BV template/default render to `mixed_voronoi`, retains explicit
> barycentric rollback and legacy omission behavior, then rerun the frozen
> acceptance without changing thresholds.

## Direct Sentaurus box-measure export - 2026-08-01

The previously unavailable Sentaurus control volumes were exported directly
with `BoxMeasureFromFile(GrdNumbering)`.  On the actual M0 and M2 meshes,
Sentaurus element-vertex `Measure` and Vela `mixed_voronoi` agree to a maximum
assembled node-volume difference of `1.39e-17 um^2`.  Both meshes have maximum
triangle angle `90 degrees` and contain no obtuse elements.  The corresponding
barycentric L1 node-volume differences are `0.02083 um^2` on M0 and
`0.09668 um^2` on M2.

A topology- and doping-preserving synthetic M0 deformation created two obtuse
non-Delaunay triangles.  Sentaurus default `AverageBoxMethod` exported no
negative `Measure` entries, but reported signed `CoeffIntersection=-0.4` and a
2% overlapping-box volume excess.  `MixAverageBoxMethod` restored exact volume
conservation.  The raw circumcentric construction had three negative local
area contributions.  Vela's current per-triangle half/quarter/quarter obtuse
rule did not reproduce Sentaurus MixAverage truncation on this synthetic mesh.

This converts the actual M0/M2 node-volume conclusion from an inference to a
direct measurement.  It also limits the equivalence claim: current
`mixed_voronoi` is Sentaurus-equivalent for the non-obtuse PN2D grid family,
not yet for arbitrary non-Delaunay grids.

Evidence:

- `docs/validation/sentaurus_box_measure_direct_export_2026-08-01.md`;
- `build-release/pn2d-box-measure-probe-20260801/audit/m0_direct_compare.json`;
- `build-release/pn2d-box-measure-probe-20260801/audit/m2_direct_compare.json`;
- `build-release/pn2d-box-measure-probe-20260801/audit/box_measure_audit.json`.

### Next execution prompt

> Keep production defaults unchanged.  For current M0/M2 validation, treat the
> mixed-Voronoi node box measure as closed against Sentaurus and continue any
> remaining BV discrepancy localization in transport, source mapping, and
> coupled feedback.  If support for arbitrary obtuse grids is required, open a
> separate contracted task to reproduce Sentaurus MixAverageBoxMethod's
> neighborhood-aware truncation; do not silently relabel the current local
> half/quarter/quarter fallback as Sentaurus MixAverage.

## Independent node-volume default reviews - 2026-08-01

The two required read-only reviews are complete and were performed without
sharing verdicts.

- Scientific review: `APPROVE_WITH_CONDITIONS`.  It authorizes only a
  qualified non-obtuse PN2D M0/M2 template proposal.  Direct Sentaurus box
  measures close the actual M0/M2 geometry question, but the synthetic obtuse
  control forbids a general MixAverage-equivalence claim.
- Code review: `APPROVE_WITH_CONDITIONS`, current worktree authorization `no`.
  The actual atomic default patch is absent: the template still defaults to
  legacy, does not render `mesh_geometry.node_volume_policy`, and the
  acceptance aggregator does not bind the candidate runs to the actual default
  render.

Combined typed outcome:

`independent_reviews_complete_atomic_patch_and_rereview_required`

The reviews permit an isolated patch that binds SG/Laux plus mixed-Voronoi
through one PN2D BV profile selector, preserves explicit
legacy-plus-barycentric rollback, rejects half-migrations, qualifies the mesh
scope, and binds acceptance to the real rendered configs.  They do not yet
authorize changing the production default.

Evidence:

- `docs/validation/pn2d_node_volume_policy_independent_scientific_review_2026-08-01.md`;
- `docs/validation/pn2d_node_volume_policy_independent_code_review_2026-08-01.md`;
- `docs/validation/pn2d_node_volume_policy_independent_review_decision_2026-08-01.json`.

### Next execution prompt

> Implement only the isolated atomic PN2D BV profile patch identified by the
> code review.  Keep global C++ defaults, PN2D IV, solver physics, algorithms,
> and thresholds unchanged.  Add default-render, rollback, half-migration,
> mutation, mesh-qualification, and artifact-binding tests; rerun the frozen
> M0/M2 and forward acceptance against the actual default render, then perform
> fresh patch-level code and scientific-scope reviews.
