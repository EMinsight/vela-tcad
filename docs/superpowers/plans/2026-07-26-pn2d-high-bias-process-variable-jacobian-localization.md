# PN2D high-bias process-variable and source-Jacobian localization plan

Date: 2026-07-26
Last updated: 2026-07-27

Status: stopped at the Task 12 nonlinear gate after completing the independent
Task 14 evidence path. Task 12 is `nonsmooth_branch_derivative`; Tasks 13 and
14 are complete. Tasks 15-19 and every self-consistent or fine PN2D run remain
unauthorized. No production-default change is authorized.

## Relationship to completed work

This plan starts from the frozen outcomes in:

- `docs/validation/pn2d_general_tri3_element_edge_avalanche_2026-07-26.md`;
- `docs/validation/pn2d_imported_state_qfp_update_2026-07-26.md`;
- `docs/validation/pn2d_tasks6_10_review_response_2026-07-26.md`;
- `docs/validation/pn2d_high_bias_process_variable_contract_2026-07-26.md`;
- `docs/validation/pn2d_exact_high_bias_oracle_2026-07-26.md`;
- `docs/validation/pn2d_element_edge_local_forward_ad_2026-07-27.md`; and
- `docs/validation/pn2d_task12_task14_closure_2026-07-27.md`.

Do not rerun the stopped Tasks 8-9 by relabeling them. The successor work must
first close the failed derivative and causality gates.

Frozen outcomes are:

- imported-state densities and P1 field/QFP-gradient vectors pass;
- native general-element alpha and documented native directed-edge carrier
  current remain insufficient observations;
- matching-support current medians pass, while P95, sign, terminal-current,
  and KCL gates fail;
- coarse source gates pass, while acute-skewed source gates fail;
- the opt-in constrained-obtuse area-conservation defect is fixed;
- before local forward AD, the source-only analytic versus central-FD Jacobian
  true-relative maximum was `0.9640948767506723`, above the frozen `1e-8`
  gate;
- the imported-state first-update outcome is
  `operator_improvement_without_qfp_causality`; and
- the production decision is `keep_opt_in_diagnostic_no_default_change`.

Execution evidence added after the Tasks 6-10 review is:

- fixed nine-slot Tri3 local forward AD evaluates and scatters all
  `3 vertices x (psi, electron QFP, hole QFP)` derivatives;
- the direct shared-source three-step FD test reaches `2.607230050e-9`, and
  unscaled/scaled scatter errors are approximately `3.4e-16`/`3.3e-16`;
- direct shared-residual and cell-contribution FD probes remove full-block
  subtraction; coarse7x3 `-10 V` closes the near-zero absolute gate, while the
  `-20 V` best nonzero true-relative result is `1.839007985e-6`;
- eight exact-flat hole-QFP edges and a minimum nonzero edge difference of
  `4.227729278e-13 V` force the real-state probe between branch crossing and
  double-precision noise, giving `nonsmooth_branch_derivative`;
- residual, AD, and cell-FD configuration identity closes at SHA-256
  `3ff38a69c6e14b011668bcbd5984ba2f7034d276338ff34b1d69f98bc893b2bd`;
- Task 13 completed as `native_process_variable_contract_available`; and
- Task 14 completed as `exact_high_bias_oracle_available`: same-cell evidence
  identifies density as the first material stage at interior hotspot cell 13.
## Scientific claim boundary

Keep these claims separate:

1. `source_jacobian_parity`;
2. `fixed_state_process_chain_localized`;
3. `first_qfp_update_causality`;
4. `coarse_self_consistent_parity`; and
5. `physical_bv_parity`.

Passing an earlier claim does not imply a later one. A fixed-state source
improvement cannot authorize a BV run without the derivative and first-update
gates.

## Priority order

| Priority | First unresolved dependency | Tasks | Authorization consequence |
|---|---|---|---|
| P0 | assembled real-state source-block parity, dependency closure, and residual/config identity | 11-12 | blocks every nonlinear claim |
| complete | same-support native Sentaurus process ordering near the knee | Task 14 | density is the first material stage on fixed interior hotspot cell 13 |
| P1 | shared Real/AD/diagnostic formula ownership and capability guards | Task 12 follow-up | prevents formula drift and premature nonlinear use |
| P1 | current/mobility/alpha support and source mapping | 15 | blocked until Tasks 12 and 14 close |
| P1 | same-residual first QFP update on both topologies | 16 | blocks continuation and sweeps |
| P2 | contact/interior and continuation branch selection | 17 | blocks self-consistent attribution |
| P2 | deterministic coarse self-consistent behavior | 18 | blocks physical fine-mesh BV |
| P3 | fine PN2D BV knee and production decision | 19 | only all prior gates may authorize review |

## Updated implementation workstreams

### P0 - freeze the canonical opt-in BV candidate

Use one validation profile:

- `generation=current_density`;
- `current_approximation=element_edge_sg_gss_laux`;
- `source_mapping_mode=element_vertex_box_measure`;
- coupled Newton; and
- the local-forward-AD source Jacobian.

This is a canonical candidate profile, not a production-default change. Do
not describe it as a unified production BV path until the Task 12 assembled
block and Task 16 cross-topology causality gates pass.

### P0 - assembled source-block Jacobian stop

The direct shared residual, local-forward-AD scatter, cell-local independent
FD contributions, scaling, constrained rows, and configuration identity have
been audited. The `-10 V` near-zero gate closes. The `-20 V` nonzero gate is
limited by exact/near-flat hole-QFP `abs` branches and returns
`nonsmooth_branch_derivative`; it is not authorized for tolerance relaxation.

### P1 - consolidate formulas and enforce capability boundaries

After the frozen P0 probes are reproducible:

1. move duplicated Real/AD Van Overstraeten, mobility, SG stabilization,
   clamp, floor, and active-branch algebra behind scalar-templated helpers;
2. distinguish `local_source_ad_validated`, `global_scatter_validated`,
   `full_block_jacobian_validated`, `qfp_causality_validated`, and
   `self_consistent_bv_authorized`;
3. reject or warn on incomplete node-local Jacobian modes in coupled Newton/BV;
4. defer ordinary SG edge-source stencil AD until the canonical element-edge
   full-block gate closes.

## Frozen constraints

1. Keep `element_edge_sg_gss_laux` opt-in.
2. Do not change a production default in Tasks 11-19. If every gate passes,
   Task 19 may propose a separate default-change review, not make the change.
3. Do not fit field, mobility, alpha, current, source, geometry, or voltage
   scales.
4. Keep Van Overstraeten and the global QFP-gradient default unchanged unless
   new direct and independent evidence falsifies them.
5. Never call a vertex-state or box-operator reconstruction a native edge
   current. Use `box_operator_reconstruction` or `operator_replay`.
6. Do not use Minimal6 as a physical BV case.
7. Use exact Sentaurus states. Do not interpolate an imported state into a
   fixed-state gate.
8. Preserve untracked `tmp/`; do not read, modify, stage, or delete it.
9. Do not commit generated `build-release/` output.
10. Use MSYS2 UCRT64. At every commit boundary run focused tests,
    `ascii_sources`, `git diff --check`, and a scoped diff.

## Bias and topology ladder

Use:

- `minimal6` only as a frozen regression and derivative micro-case;
- `coarse7x3` as the high-bias knee, contact/interior, and controlled
  self-consistent case;
- `skewed_tri3` and `constrained_obtuse` as geometry diagnostics; and
- `fine_pn2d` only after Task 18 passes.

The initial exact high-bias lattice is:

`-10.0, -18.0, -19.0, -19.5, -19.8, -19.9, -19.95, -20.0 V`.

Any finer lattice must be declared and hashed before comparing Vela and
Sentaurus. Preserve every accepted state and the first rejected transition.

## Task 11 - freeze the successor schema and RED contracts

### Purpose

Prevent a diagnostic definition from changing after the first failed
derivative or process-chain stage is observed.

### Actions

1. Add a versioned schema for:
   - derivative contributions;
   - native versus reconstructed observation labels;
   - process-variable stage values;
   - topology/contact class;
   - residual and Jacobian configuration hashes;
   - exact state/deck/TDR/mesh/parameter/release provenance; and
   - finite, zero, below-floor, nonsmooth, unsupported, or invalid status.
2. Freeze central-FD step-selection and convergence checks for
   `psi`, electron QFP, and hole QFP.
3. Freeze:
   - nonzero Jacobian true-relative gate `<= 1e-8`;
   - near-zero absolute gate `<= 1e-12`;
   - source/residual identity gate `<= 1e-12`; and
   - derivative contribution-sum closure `<= 1e-12`.
4. Add RED tests for:
   - a missing dependency contribution;
   - different residual/Jacobian configurations;
   - a reconstructed current labeled native;
   - an exact zero converted into a finite dex value;
   - a denominator using `max(1, norm)` instead of true relative error; and
   - a non-exact or hash-mismatched state.

### Exit gate

- RED fails for the intended missing schema or contract.
- Failure is not caused by locale, newline, path, or absent old build output.

### Typed outcomes

- `successor_red_contract_frozen`; or
- `schema_contract_red_failed`.

### Commit boundary

Commit only the schema, verifier contract, and RED tests.

## Task 12 - isolate and close the source-only Jacobian

### Purpose

Identify the first analytic dependency that does not match the independently
evaluated central finite difference.

### Current status

Stopped with typed outcome `nonsmooth_branch_derivative`. The local forward-AD
source, global scatter, direct shared residual, cell-local FD contributions,
scaling, constrained rows, and configuration identity are validated. The
historical `0.9640948767506723` and `1.060161811e-2` results are pre-direct
probe baselines. The final `-20 V` best true-relative result is
`1.839007985e-6`, above `1e-8`; smaller same-branch steps enter roundoff.
### Entry condition

Task 11 RED is accepted.

### Actions

1. Write the source chain explicitly for both carriers:

   `state -> density -> low-field mobility -> high-field mobility/drive ->`
   `SG edge current -> cell current vector -> abs(J) -> alpha(drive) ->`
   `element source -> node residual`.

2. Record analytic and central-FD contributions to every affected
   `psi`, electron-QFP, and hole-QFP block for:
   - carrier-statistics and intrinsic-density terms;
   - SG Bernoulli and potential/QFP differences;
   - low-field mobility;
   - high-field mobility versus QFP-gradient or electric-field drive;
   - current-vector reconstruction;
   - `abs(J)`, carrier sign, zero floor, and active branch;
   - Van Overstraeten `alpha(F)` and its driving-force derivative;
   - element-edge source distribution;
   - element-to-node accumulation; and
   - continuity scaling.
3. Require residual assembly, diagnostic replay, and analytic Jacobian to
   emit identical mobility, driver, source-map, clamp, and floor hashes.
4. Preserve the RED case that reproduced the historical
   `0.9640948767506723` failure, and use the frozen coarse7x3 `-20 V`
   assembled real-state probe as the current nonzero gate.
5. If one dependency is missing or inconsistent, make the smallest opt-in
   patch and rerun the same RED case before adding broader cases.
6. Test FD convergence with at least three symmetric step sizes. A derivative
   crossing a sign, floor, or truncation branch is classified as nonsmooth;
   it is not made GREEN by relaxing the tolerance.
7. Cover acute, constrained-obtuse, reversed-orientation, contact-adjacent,
   and interior cases for both carriers.
8. Evaluate the shared source residual directly; do not infer the small source
   block by subtracting numerically dominant full-continuity blocks.
9. Require signed edge flux, current vector, alpha, and source diagnostics to
   come from, or be checked against, the same shared source core.

### Exit gate

- Nonzero analytic/central-FD true-relative maximum: `<= 1e-8`.
- Near-zero analytic/central-FD absolute maximum: `<= 1e-12`.
- Analytic contribution sum closes to the assembled analytic block:
  `<= 1e-12`.
- Independent FD contribution sum closes to the full source-only FD block:
  `<= 1e-12`.
- Diagnostic source and residual source are identical: `<= 1e-12`.
- Residual and Jacobian configuration hashes are identical.
- No unrelated Poisson, SRH, statistics, mobility, or default formula changes.

### Typed outcomes and stop condition

Allowed outcomes:

- `source_jacobian_dependency_identified_and_closed`;
- `incomplete_analytic_derivative`;
- `nonsmooth_branch_derivative`;
- `residual_jacobian_configuration_mismatch`; or
- `jacobian_gate_failed`.

Stop all nonlinear authorization if the first outcome is not achieved.
Preserve the dependency-level evidence and do not enter Tasks 16-19.

### Commit boundary

Commit focused tests first where practical, then only the minimal opt-in
derivative/configuration patch.

## Task 13 - prove the Sentaurus process-variable export contract

### Purpose

Determine which high-bias quantities are directly observable before using
them to explain the `-19` to `-20 V` current turn.

### Current status

Complete: `native_process_variable_contract_available`. The accepted exact
coarse7x3 state is `-19.95 V` on Sentaurus Device `O-2018.06-SP2`. This
completion authorizes Task 14 evidence only and does not override Task 12.

### Entry condition

Task 11 is complete. This evidence task may proceed if Task 12 stops, but it
cannot authorize a Vela nonlinear candidate by itself.

### Actions

1. On `coarse7x3`, make a one-state exact export probe at `-19.95 V`; repeat
   on `fine_pn2d` only after the coarse schema is proven.
2. Request ordinary Plot/TDR fields already supported by the deck:
   - potential and electron/hole QFP;
   - electron/hole density;
   - electric field and element vector;
   - electron/hole QFP-gradient element vectors;
   - electron/hole mobility and velocity;
   - electron/hole/total current-density vectors;
   - electron/hole avalanche alpha;
   - carrier-split avalanche generation;
   - electron/hole and mean ionization integrals;
   - doping, space charge, and SRH recombination.
3. Request CurrentPlot observations:
   - terminal currents;
   - carrier-split source integrals;
   - element-local `ReadCoefficient`; and
   - element-local-vertex `ReadMeasure`.
4. Verify field centering, component count, units, region/cell/node mapping,
   carrier sign, and exact bias. In particular, require two current-vector
   components and scalar electron/hole alpha on the accepted TDR.
5. Repeat the raw export into two independent roots and compare SHA-256
   ledgers.
6. Probe `CNormPrint` or `NewtonPlot` only if the installed release has a
   locally documented minimal syntax. Do not claim full residual or Jacobian
   export from ordinary Plot/CurrentPlot.
7. Preserve negative evidence:
   - `/Edge` Plot and ReadFlux rejection do not become native edge data;
   - native directed-edge carrier current remains unsupported unless directly
     exported; and
   - operator replay remains reconstructed.

### Exit gate

- Every accepted field has an explicit native/reconstructed/unsupported label.
- Exact bias, release, deck, mesh, TDR, and parameter hashes are sealed.
- Current vectors have the required component count.
- Electron/hole alpha and requested process variables are present or each has
  a typed unsupported result.
- Two raw roots agree byte-for-byte after normalized path metadata.

### Typed outcomes

- `native_process_variable_contract_available`;
- `partial_native_process_variable_contract`;
- `main_mesh_process_variable_contract_unavailable`; or
- `insufficient_native_observation`.

Only the first two outcomes may enter Task 14, and the second must carry an
explicit list of claims it cannot support.

## Task 14 - generate the exact high-bias Sentaurus oracle

### Purpose

Resolve where the Sentaurus current begins to accelerate and which native
process variable changes first.

### Current status

Complete: `exact_high_bias_oracle_available`. All
`9 branches x 8 exact biases x 2 roots = 144` observations and paired
provenance are accepted. Same-cell contact/interior and active-region analysis
fixes the `-20 V` generation hotspot at interior cell 13 and identifies density
as the first material stage from `-19 V` to `-20 V`.
### Entry condition

Task 13 has an accepted native process-variable contract.

### Actions

1. Generate every exact state in the frozen high-bias lattice on two
   independent raw roots.
2. Use otherwise identical branches:
   - production implicit avalanche;
   - explicit `GradQuasiFermi`;
   - explicit electric-field drive;
   - avalanche disabled;
   - `AvalDensGradQF`;
   - paired high-field-saturation-disabled controls; and
   - contact-selector control only where its sole difference is proven.
3. Preserve the exact terminal current curve and compute adjacent-voltage
   growth without interpolating states.
4. At every state, export the Task 13 process variables and compute
   contact-adjacent/interior and active-region summaries.
5. Detect the first material change in dependency order:

   `psi/QFP -> density -> field/QFP gradient -> mobility/velocity ->`
   `current -> alpha/ion integral -> source -> terminal current`.

6. Seal solver log, accepted/rejected transition, deck, release, TDR, mesh,
   model, and Math-option provenance.
7. Reuse the accepted matrix to align density, mobility, current, alpha,
   generation, and terminal response on the same cells or hotspot; do not use
   maxima from different locations to claim process ordering.

### Exit gate

- Every lattice row is exact, converged, and hash-paired.
- Paired branches differ only by the declared selector or model control.
- The first high-bias departure is independently reproduced from raw fields.
- No failed or interpolated state enters a gate.

### Typed outcomes

- `exact_high_bias_oracle_available`;
- `partial_high_bias_oracle`;
- `sentaurus_branch_provenance_mismatch`; or
- `high_bias_oracle_unavailable`.

Only the first outcome may enter Task 15.

### Stop condition

Stop if a branch changes undeclared physics, a high-field-saturation control
is unpaired, or current/alpha component schemas change across bias.

## Task 15 - localize current, mobility, alpha, and source support

### Purpose

Find the first non-closing fixed-state stage without conflating a Sentaurus
box current with a Vela SG edge current.

### Entry condition

Task 12 returns `source_jacobian_dependency_identified_and_closed` and Task 14
returns `exact_high_bias_oracle_available`. Task 14 now satisfies its side of
the condition, but the current Task 12 `nonsmooth_branch_derivative` outcome
does not. Task 15 is therefore not entered.

### Actions

For each exact imported state, compare in stages:

1. imported `psi/QFP -> n/p`;
2. P1 electric field and QFP gradients;
3. low-field and high-field mobility, including field support;
4. element-local directed SG edge-current `operator_replay`;
5. GSS/Laux cell-current-vector reconstruction;
6. native Sentaurus element current-density vector;
7. alpha and ionization-integral drive;
8. element-vertex source using `ReadMeasure`;
9. physical-node source accumulation;
10. device source integral and terminal current.

Run staged substitutions of one stage at a time. Report separately:

- electron and hole;
- contact-adjacent and interior cells;
- acute, right, and obtuse geometry;
- all-three-edge sensitivity;
- magnitude, component, angle, and sign errors;
- active-range and near-zero tails; and
- native, reconstructed, and unsupported observations.

Use terminal-current and KCL gates only on legitimate global fields and
terminal observations. Do not infer native edge-current parity from a local
operator replay.

### Exit gate

- Previously passed density and P1 gates remain unchanged.
- Constant-field vector and geometric support closure: `<= 1e-12`.
- Source identity: `<= 1e-12`.
- `ReadMeasure` versus CurrentPlot integral: `<= 1e-10` relative.
- Full reconstructed Sentaurus source maximum: `<= 5e-3 dex`.
- Matching-support current median/P95: `<= 0.05/0.15 dex`.
- Nonzero carrier sign agreement: 100%.
- Vela fixed-state source median/maximum: `<= 0.10/0.30 dex` in the
  avalanche-active range.
- Legitimate terminal-current closure: `<= 2e-7` relative.
- Legitimate internal total-current KCL: `<= 1e-8` relative.
- One earliest failed stage is identified for every failed topology/carrier
  class.

### Typed outcomes

- `fixed_state_process_chain_closed`;
- `mobility_support_difference`;
- `element_current_support_difference`;
- `source_mapping_geometry_difference`;
- `contact_interior_model_difference`;
- `proprietary_operator_difference`; or
- `insufficient_native_observation`.

Only `fixed_state_process_chain_closed`, or a directly proven minimal Vela
defect that is patched and then closes the same gates, may authorize Task 16.

## Task 16 - re-establish first-QFP-update causality

### Purpose

Test whether the authorized source/current correction improves the actual
continuity residual and first QFP update on both diagnostic topologies.

### Entry condition

Task 12 Jacobian gates pass and Task 15 authorizes one opt-in candidate.

### Actions

1. At exact `-10`, `-19`, `-19.95`, and `-20 V` imported states, assemble
   baseline and candidate electron/hole continuity terms on Minimal6 and
   coarse7x3.
2. Record SG divergence, SRH, carrier-split avalanche source, boundary term,
   physical and normalized residual, row scale, and final residual.
3. Compute carrier-only and coupled first Newton updates with identical
   configuration, damping, clamps, and linear solver.
4. Record the first bias/node/carrier/term that changes the update and compare
   it with the Sentaurus-to-Vela QFP error direction.
5. Keep avalanche-off and SRH-off paired controls for production and candidate
   configurations.

### Exit gate

- Residual decomposition closure: `<= 1e-12`.
- Boundary rows remain identical to floating-point precision.
- Task 12 Jacobian gates remain satisfied.
- First updates exist for every required bias, topology, and carrier.
- The same candidate improves the same residual and first-update QFP metric on
  both topologies without worsening the other carrier or a lower bias.
- Independent verification derives the authorization outcome rather than
  reading a preselected label.

### Typed outcomes and stop condition

- `source_support_causes_qfp_update`;
- `current_or_mobility_causes_qfp_update`;
- `boundary_or_contact_model_difference`;
- `srh_or_source_term_difference`;
- `operator_improvement_without_qfp_causality`; or
- `proprietary_model_difference`.

Only one of the first two outcomes with all exit gates may enter Task 17.
Otherwise stop before any self-consistent sweep.

## Task 17 - localize continuation and branch selection

### Purpose

Determine whether the remaining high-bias difference is a nonlinear path
effect after the fixed-state operator and first update have been authorized.

### Entry condition

Task 16 passes on both topologies.

### Actions

1. Hold physics, mesh, contacts, scaling, tolerances, and candidate operator
   hashes fixed.
2. Compare:
   - cold versus exact warm restart;
   - at least two predeclared bias step schedules;
   - carrier-only versus coupled diagnostic first steps; and
   - contact-adjacent versus interior first-failure nodes.
3. Preserve every accepted checkpoint, rejected transition, Newton residual,
   damping factor, clamp activation, and linear-solver status.
4. Repeat every branch twice and require deterministic state/output hashes.
5. Treat continuation as a diagnosis, not a fitted voltage or current
   correction.

### Exit gate

- Identical physics/configuration hashes across schedules.
- Duplicate runs are deterministic.
- A schedule-independent fixed point agrees within the frozen Task 18 state
  targets, or the first branch-dependent transition is localized.
- No schedule is selected solely because it hides an operator mismatch.

### Typed outcomes

- `continuation_invariant_candidate`;
- `continuation_branch_difference`;
- `contact_adjacent_first_failure`;
- `interior_first_failure`; or
- `solver_first_failure`.

Only `continuation_invariant_candidate` may enter Task 18.

## Task 18 - deterministic coarse self-consistent sweep

### Purpose

Test whether the opt-in candidate preserves the coarse high-bias branch and
reproduces the Sentaurus current turn without regression.

### Entry condition

Tasks 12, 15, 16, and 17 authorize the same unchanged candidate.

### Actions

1. Run baseline and candidate twice on coarse7x3 using the frozen exact bias
   lattice and the authorized continuation schedule.
2. Compare in dependency order:

   `psi -> QFP -> n/p -> mobility -> current -> alpha -> source ->`
   `terminal current`.

3. Re-run Tasks 12 and 15 at every accepted self-consistent checkpoint.
4. Report the `-19` to `-20 V` interval explicitly:
   - current at each exact state;
   - adjacent-state growth ratio;
   - first `1.5x` and `2.0x` one-volt-equivalent markers;
   - first variable departing before each marker; and
   - the first rejected transition.

### Exit gate

| Quantity | Target |
|---|---:|
| electrostatic potential maximum | `1e-6 V` |
| electron/hole QFP median/P95 | `0.01/0.025 V` |
| electron/hole density median/P95 | `0.10/0.25 dex` |
| matched-support mobility median/P95 | `0.05/0.20 dex` |
| matching-support operator-replay current median/P95 | `0.10/0.25 dex`, 100% nonzero sign; not a native-edge claim |
| total terminal current median | `0.10 dex` |
| active-range impact-source median | `0.30 dex` |
| internal total-current KCL | `1e-8` relative |

Additionally:

- no accepted low-bias state regresses relative to baseline;
- duplicate hashes and first rejected transition agree;
- source/Jacobian and terminal/KCL closures remain valid; and
- the current-turn conclusion is based on exact rows, not interpolation.

### Typed outcomes and stop condition

- `coarse_self_consistent_parity_passed`;
- `operator_improves_but_qfp_misses`;
- `fixed_state_only_improvement`;
- `continuation_branch_difference`;
- `solver_first_failure`; or
- `model_difference`.

Only the first outcome may enter Task 19.

## Task 19 - fine PN2D BV, decision, validation, and review

### Purpose

Evaluate physical BV only after the candidate survives the coarse nonlinear
authorization chain.

### Entry condition

Task 18 returns `coarse_self_consistent_parity_passed`.

### Actions

1. Prove the Task 13 native export schema on `fine_pn2d`.
2. Run baseline and the unchanged authorized candidate with identical mesh,
   contacts, parameters, bias stepping, continuation, tolerances, and stop
   conditions.
3. Preserve exact rows and the first rejected transition.
4. Compare:
   - `1.5x` and `2.0x` current-growth knee markers;
   - maximum log-current error over `-20..-10 V`;
   - QFP, mobility, field, alpha, ionization integral, source, terminal
     current, and KCL near each marker; and
   - `latent_turning_point` versus `physics_magnitude_gap`.
5. Report a missing marker as `undefined` and type it
   `knee_not_observed`.
6. Run focused tests, general-mesh regressions, two-root independent
   verification, full Release build/CTest, `ascii_sources`, manifest/schema
   validation, `git diff --check`, and scoped diffs.
7. Request independent scientific review and code review.

### Exit gate

- Both Sentaurus markers exist.
- Candidate `1.5x` and `2.0x` marker biases are each within `1.0 V`.
- Maximum knee-window log-current error: `<= 0.30 dex`.
- No unclassified failure occurs before the Sentaurus knee.
- Terminal-current and KCL closures remain valid.
- Ionization-integral, source, and current changes have one consistent causal
  direction.

### Typed outcomes

- `physical_bv_parity_passed`;
- `knee_not_observed`;
- `latent_turning_point`;
- `physics_magnitude_gap`;
- `solver_first_failure`; or
- `model_difference`.

Only `physical_bv_parity_passed` may propose a separate default-change
review.

### Decision

| Evidence | Allowed decision |
|---|---|
| Any prior gate fails | keep opt-in diagnostic; preserve typed outcome |
| Coarse passes but fine BV fails | keep opt-in; classify physical BV gap |
| Every gate and both reviews pass | propose a separate default-change review |
| Any result needs a fitted scale | `model_difference`; no production patch |

### Commit structure

Keep commits scoped:

1. Task 11 schema and RED contracts;
2. Task 12 focused Jacobian tests and minimal opt-in patch, if authorized
   (local forward AD is committed; assembled-block closure remains);
3. Tasks 13-15 export/analyzer/verifier/report (Task 13 is committed, Task 14
   is partial, and Task 15 is not entered);
4. Task 16 first-update causality evidence;
5. Task 17 continuation evidence;
6. Task 18 coarse self-consistent evidence;
7. Task 19 fine BV evidence and review responses.

Never stage unrelated worktree changes or generated simulation roots.

## Global stop conditions

Stop and preserve partial evidence if:

1. residual and Jacobian use different mobility, driver, source mapping,
   clamp, or floor;
2. a local source or scatter test passes while the frozen assembled real-state
   source block fails;
3. a derivative crosses a nonsmooth branch without being typed;
4. paired Sentaurus branches differ in undeclared physics or provenance;
5. topology, cell, component, unit, carrier sign, or edge orientation is
   ambiguous;
6. a reconstructed quantity is labeled native;
7. source identity or `ReadMeasure` integration does not close;
8. a zero or low-signal tail is converted into a finite active-range dex
   claim;
9. more than one causal factor changes in an isolated comparison;
10. a fitted scale or proprietary inferred value enters a candidate;
11. a first-update result from one topology is generalized to another;
12. continuation is tuned before fixed-state and first-update gates pass;
13. a coarse failure is bypassed with a fine run;
14. Minimal6 is called a physical BV case; or
15. a default change is made inside this plan.

## Final success definition

This plan succeeds with either:

1. an independently reviewed chain from source-Jacobian parity through
   native process variables, fixed-state support, first QFP update, coarse
   self-consistency, and physical BV; or
2. a bounded typed outcome naming the first non-closing dependency and
   preserving the current accurate formulas and opt-in/default boundary.
