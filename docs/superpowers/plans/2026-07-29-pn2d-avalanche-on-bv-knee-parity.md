# PN2D avalanche-on BV curve and knee parity plan

Date: 2026-07-29

Status: engineering WP0-WP7 executed; Tasks 5-6 returned `insufficient_observation`, so Tasks 7-11 are not authorized.

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
| Tasks 5-6 | WP3-WP5 records and WP7 paired process analyzer | `insufficient_observation`: matching Vela exact-lattice process manifest absent |
| Tasks 7-11 | causal authorization from Tasks 4-6 | prohibited; no two-bias causal stage |

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
| M0 passes but dose-preserving refinement fails | keep opt-in; classify mesh dependence |
| Curves match only after scale/shift/fitting | reject; classify `model_difference` |
| All tasks and both reviews pass | propose a separate production-default change |

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

> Execute
> `docs/superpowers/plans/2026-07-29-pn2d-avalanche-on-bv-knee-parity.md`
> from WP2.1/Task 1's incomplete exact-lattice gate, then implement Task 2.
> Preserve the sealed deterministic Vela failure at
> `-19.692187499999644 -> -19.693749999999643 V` as an explicit baseline
> result. Strictly obey every entry condition, acceptance criterion, decision
> gate, and stop condition. Keep BV on `net_doping`, forward IV on
> `cell_reconstructed_total_impurity`, avalanche drive on
> `quasi_fermi_gradient`, and Van Overstraeten coefficients unchanged. Report
> machine-readable evidence after every decision gate before proceeding. Do
> not stage `tmp/` or commit generated simulation outputs.
