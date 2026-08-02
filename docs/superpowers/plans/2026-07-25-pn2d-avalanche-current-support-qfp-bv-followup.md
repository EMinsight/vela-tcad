# PN2D avalanche current-support, QFP, and BV follow-up plan

Date: 2026-07-25

Status: executed through the mandated Task 6/7 stop conditions on 2026-07-26;
Tasks 8-9 were not entered, and no production-default change is authorized.

Execution ledger: `docs/validation/pn2d_tasks6_10_review_response_2026-07-26.md`.

Successor plan:
`docs/superpowers/plans/2026-07-26-pn2d-high-bias-process-variable-jacobian-localization.md`.

## Relationship to completed work

This plan begins after completion of the Minimal6 box-edge current-factor,
element-edge fixed-state, and self-consistent audits, together with:

- the Sentaurus avalanche-default correction; and
- the final independent scientific and code reviews.

It supersedes only the former immediate-next-action statement. All completed
Minimal6 evidence remains frozen and must be rerun as regression evidence, not
reinterpreted or refitted.

The existing uncommitted naming-only edits in
`include/vela/physics/ImpactIonizationModel.h` and
`src/physics/ImpactIonizationModel.cpp` are outside this plan. Preserve them
or commit them separately before starting Task 1. Do not stage `tmp/`.

## Objective

Determine whether Vela can reproduce the Sentaurus avalanche source on
general Tri3 meshes by aligning:

1. element-local Scharfetter-Gummel carrier currents;
2. GSS/Laux current-vector reconstruction;
3. element-vertex source measures;
4. contact versus interior-element avalanche driving-force selection;
5. residual and Jacobian assembly; and
6. the self-consistent QFP branch and physical BV knee.

The plan must distinguish three claims:

- `fixed_state_operator_parity`;
- `self_consistent_state_parity`; and
- `physical_bv_parity`.

Passing an earlier claim does not imply a later one.

## Frozen evidence

The following observations are inputs and may not be refitted:

- Minimal6 exact lattice:
  `mirror/sketch x -1..-20 V`, 40 states.
- Imported Sentaurus `psi/phin/phip` plus Vela BGN reproduces all node
  densities within `4.42618e-6 dex`.
- Vela P1 `-grad(psi)` matches native Sentaurus element electric-field vectors
  to `2.959312e-16` maximum relative residual.
- Minimal6 effective avalanche drive is electric field because every triangle
  touches a contact. The global Sentaurus drift-diffusion default remains
  `GradQuasiFermi`.
- Electric-field-driven Van Overstraeten coefficient errors are at most
  `1.55e-9 dex` for electrons and `2.56e-9 dex` for holes.
- Minimal6 element-edge fixed-state source-integral median errors are
  `0.0548177 dex` electron and `0.0517363 dex` hole.
- A Sentaurus-native box-operator replay closes the source integral with
  `4.73942e-5 dex` median and `4.29237e-4 dex` maximum error.
- The opt-in self-consistent element-edge source error is `1.23508 dex` over
  the full low-signal lattice, `0.154333 dex` over `-16..-20 V`, and
  `0.0494843 dex` at `-20 V`.
- The prior triangle source proxy differs by about `10.76 dex` and is not a
  Sentaurus parity operator.
- Self-consistent electrostatic potential already passes, while electron and
  hole QFP medians remain approximately `0.0515 V` and `0.0591 V`.
- Sentaurus native directed-edge current is unavailable. Any reconstructed
  Sentaurus edge current must remain labeled
  `box_operator_reconstruction`.

## Mesh ladder

Use four mesh levels with separate scientific roles:

1. `minimal6`
   - regression-only;
   - all elements touch contacts;
   - right-triangle zero-diagonal behavior remains frozen.
2. `coarse7x3`
   - reuse `reference_tcad/pn2d_sentaurus2018_coarse7x3`;
   - discriminate contact-adjacent and interior-element behavior;
   - support controlled self-consistent and coarse BV experiments.
3. `skewed_tri3`
   - create only if the regenerated coarse7x3 mesh lacks at least one acute
     scalene and one obtuse triangle with nonzero three-edge support;
   - diagnostic-only Sentaurus/Vela operator mesh.
4. `fine_pn2d`
   - reuse `reference_tcad/pn2d_sentaurus2018`;
   - physical BV curve and knee comparison only after Tasks 1-8 pass.

Do not use Minimal6 to authorize a general-mesh production default.

## Expected work products

Likely new files:

- `scripts/run_pn2d_general_tri3_sentaurus_avalanche_controls_vm.py`;
- `scripts/diagnose_pn2d_general_tri3_element_edge_avalanche.py`;
- `scripts/verify_pn2d_general_tri3_element_edge_avalanche.py`;
- `scripts/diagnose_pn2d_imported_state_qfp_update.py`;
- `scripts/verify_pn2d_imported_state_qfp_update.py`;
- `tests/regression/test_pn2d_general_tri3_element_edge_avalanche.py`;
- `tests/regression/test_pn2d_imported_state_qfp_update.py`;
- `docs/validation/pn2d_general_tri3_element_edge_avalanche_2026-07-XX.md`;
- `docs/validation/pn2d_imported_state_qfp_update_2026-07-XX.md`; and
- deterministic evidence roots under `build-release/`.

Reuse and extend where possible:

- `tests/test_element_edge_gss_laux_avalanche.cpp`;
- `tests/test_cell_reconstructed_avalanche.cpp`;
- `tests/test_impact_ionization.cpp`;
- `scripts/diagnose_pn2d_minimal6_element_avalanche_replay.py`;
- `scripts/diagnose_pn2d_minimal6_element_edge_gss_laux_fixed_state.py`;
- `scripts/diagnose_pn2d_minimal6_phase_e_continuity_residual.py`;
- `scripts/diagnose_pn2d_bv_knee_shape.py`; and
- `scripts/diagnose_pn2d_bv_knee_gap.py`.

Do not duplicate a completed Minimal6 analyzer merely to rename its output.

## Task 1 - freeze the general-mesh schema and create RED tests

### Purpose

Prevent contact fallback, element order, edge orientation, geometric zero, or
driver selection from being silently mixed.

### Actions

1. Add a regression test for a new schema
   `pn2d_general_tri3_element_edge_avalanche/v1`.
2. Require every record to contain:
   - case, topology, bias, carrier;
   - Vela cell id and Sentaurus region-cell id;
   - node ids and local-edge orientation;
   - contact-adjacent and interior-element flags;
   - triangle angles, signed area, and orientation;
   - `ReadCoefficient` and `ReadMeasure`;
   - electric-field and QFP-gradient source hashes;
   - low-field and final mobility source hashes;
   - avalanche coefficient drive and current-density approximation;
   - native/reconstructed/unsupported observation label;
   - zero, below-floor, finite, or valid status; and
   - units for every physical column.
3. Require exact state lists. Interpolation is forbidden in fixed-state
   operator gates.
4. Add negative tests for:
   - wrong carrier sign;
   - wrong cell permutation;
   - contact fallback mislabeled as global default;
   - a zero coefficient converted to a finite dex value;
   - a reconstruction mislabeled as native;
   - mismatched TDR or `models.par` hashes; and
   - mixed Sentaurus releases.
5. Make the initial focused test fail because the general-mesh schema and
   output do not yet exist.

### RED gate

- Failure is caused by the missing schema/implementation.
- Failure is not caused by newline, locale, missing old build output, or path
  formatting.

### Commit boundary

Commit only the schema contract and RED tests.

## Task 2 - regenerate and classify the Sentaurus oracle meshes

### Purpose

Obtain an oracle containing both contact-adjacent and interior cells and
enough triangle geometry to test all edge weights.

### Actions

1. Regenerate coarse7x3 on Sentaurus O-2018.06-SP2.
2. Start with exact biases `-1`, `-10`, and `-20 V`.
3. Generate these branches with otherwise identical decks:
   - implicit default avalanche;
   - explicit `GradQuasiFermi`;
   - explicit `ElectricField`;
   - explicit `GradQuasiFermi` plus
     `ComputeGradQuasiFermiAtContacts=UseQuasiFermi`;
   - `ElectricField` plus `UseQuasiFermi`, used only as a selector control;
   - mobility-isolated electric-field versus GradQF avalanche branches with
     HighFieldSaturation disabled in both; and
   - `AvalDensGradQF` as a current-density-approximation control.
4. Export:
   - node `psi`, electron/hole QFP, density, mobility, alpha, and generation;
   - element electric field, electron/hole QFP gradients, mobility, and
     current-density vectors;
   - element-local `ReadCoefficient`;
   - element-vertex `ReadMeasure`;
   - carrier-split CurrentPlot source integrals;
   - terminal currents; and
   - runtime logs containing selected drive and Math options.
5. Compute all cell angles and contact adjacency after download.
6. If coarse7x3 lacks either:
   - an interior cell;
   - an acute scalene cell with all three positive supports; or
   - an obtuse cell exercising truncation,
   create the smallest `skewed_tri3` diagnostic deck that supplies the
   missing class, then repeat the three-bias export.
7. Produce two independent raw roots and verify remote/local SHA-256 ledgers.

### Exit gate

- Every requested branch and exact bias is present.
- At least one interior element is observed.
- Required angle classes are observed, or the skewed mesh is generated.
- Implicit/explicit driver equivalence is classified separately for contact
  and interior elements.
- Static TDR and `models.par` hashes match across paired branches.
- No failed or silently interpolated state enters the accepted lattice.

### Stop condition

Stop if the driver branches change any physics other than the declared
selector, or if HighFieldSaturation is disabled in only one mobility-isolation
branch.

## Task 3 - fixed imported-state upstream replay

### Purpose

Prove the state transformation before evaluating current or avalanche source.

### Actions

For every accepted state:

1. Import Sentaurus `psi/phin/phip`.
2. Recompute Vela `n,p` with the active Old-Slotboom/BGN statistics.
3. Recompute P1 element vectors:
   - `E = -grad(psi)`;
   - `Gn = -grad(phin)`; and
   - `Gp = -grad(phip)`.
4. Compare P1 vectors with native Sentaurus element vectors by:
   - component error;
   - magnitude error;
   - angle error; and
   - contact versus interior class.
5. Recompute low-field and high-field mobility without parameter fitting.
6. Evaluate Van Overstraeten alpha with:
   - the mesh-effective Sentaurus selector;
   - forced electric field;
   - forced QFP gradient; and
   - low-density interpolation where configured.
7. Record electron and hole alpha separately.

### Exit gate

- Electric-field vector maximum relative error: `<= 1e-12`.
- Imported-state density maximum error: `<= 1e-4 dex`.
- Matching-driver alpha maximum error: `<= 1e-6 dex` in the avalanche-active
  range.
- All sign and angle conventions are explicit.
- Electron native-QFP-gradient output differences remain typed as output or
  proprietary element-evaluation differences unless independently explained.

### Decision

Failure before current is evaluated blocks every later source-parity claim.
It does not authorize a fitted field scale.

## Task 4 - replay element-local SG current and cell vectors

### Purpose

Align the current support used by the avalanche source on general triangles.

### Actions

1. Recompute every element-local directed SG edge current from endpoint
   state, variable intrinsic density, and declared element mobility.
2. Preserve both raw edge current density and
   `ReadCoefficient * edge current`.
3. Reconstruct cell current vectors with:
   - current Vela triangle proxy;
   - GSS/Laux element-edge weighting;
   - active-edge exact control where algebraically defined;
   - Charon Whitney/HCurl control; and
   - Genius tangent least-squares control.
4. On acute cells, perturb each edge independently and require all three
   positive supports to affect the vector.
5. On obtuse cells, verify the declared truncation policy and exact area
   conservation.
6. On right triangles, retain the exact zero-hypotenuse regression.
7. Compare:
   - vector components and magnitude;
   - current versus negative-QFP-gradient angle;
   - contact carrier currents;
   - total terminal current; and
   - internal KCL.
8. Label Sentaurus edge values
   `box_operator_reconstruction`, never `native_edge_current`.

### Exit gate

- Current-vector constant-field recovery: `<= 1e-12` relative.
- Area/weight closure: `<= 1e-12` relative.
- Matching-support cell-current median error: `<= 0.05 dex`.
- Matching-support cell-current P95 error: `<= 0.15 dex`.
- Carrier signs: 100% on nonzero supports.
- Total terminal-current closure: `<= 2e-7` relative.
- Internal total-current KCL: `<= 1e-8` relative.
- Near-zero edge tails are reported in both dex and absolute current.

### Typed outcomes

- `element_edge_current_parity`;
- `mobility_limited_current_difference`;
- `support_limited_current_difference`;
- `proprietary_operator_difference`; or
- `insufficient_native_observation`.

## Task 5 - replay avalanche source support and driver selection

### Purpose

Separate alpha, current-vector, and source-mapping errors.

### Actions

1. For each element and carrier compute:

   `qG = alpha * abs(J_cell) * element_vertex_measure`.

2. Compare independently:
   - alpha;
   - cell current magnitude;
   - element-vertex source;
   - accumulated physical-node source;
   - per-cell source integral;
   - carrier-specific device integral;
   - total device integral; and
   - Sentaurus CurrentPlot integral.
3. Run staged replacements:
   - Vela alpha plus Vela current/support;
   - Sentaurus alpha plus Vela current/support;
   - Vela alpha plus reconstructed Sentaurus current/support;
   - full reconstructed Sentaurus operator; and
   - each driver-control branch.
4. Report contact-adjacent and interior cells separately.
5. Preserve low-signal states as typed diagnostics but exclude them from
   avalanche-active dex gates.
6. Confirm `AvalDensGradQF` changes current-density approximation rather than
   the alpha driving-force selector.

### Exit gate

- Source identity maximum relative error: `<= 1e-12`.
- ReadMeasure versus CurrentPlot integral: `<= 1e-10` relative.
- Full reconstructed Sentaurus operator source maximum error:
  `<= 5e-3 dex`.
- Vela fixed-state element-edge source median error:
  `<= 0.10 dex` in the active range.
- Vela fixed-state element-edge source maximum error:
  `<= 0.30 dex` in the active range.
- Contact and interior driver behavior agrees with the documented Sentaurus
  selector contract.

### Decision

If alpha closes and source changes only when current support is replaced, keep
the Van Overstraeten formula unchanged.

## Task 6 - generalize the C++ operator and seal residual/Jacobian parity

### Entry condition

Tasks 3-5 must identify a reproducible Vela implementation or support defect.
If the existing opt-in implementation already passes, add only tests and
evidence; do not rewrite it.

### Actions

1. Start with RED tests in:
   - `tests/test_element_edge_gss_laux_avalanche.cpp`;
   - `tests/test_cell_reconstructed_avalanche.cpp`; and
   - `tests/test_impact_ionization.cpp`.
2. Cover:
   - acute scalene cell;
   - obtuse cell;
   - forward/reverse orientation;
   - contact and interior driver selection;
   - electron/hole signs;
   - all-three-edge sensitivity;
   - zero-hypotenuse regression;
   - source identity;
   - residual/Jacobian use of identical mobility and driver configuration;
   - carrier-specific alpha; and
   - near-zero source behavior.
3. If required, make the smallest patch in:
   - `include/vela/equation/AssemblerUtils.h`;
   - `src/equation/CoupledDDAssembler.cpp`; or
   - their directly related headers.
4. Keep `element_edge_sg_gss_laux` opt-in.
5. Do not change the production default in this task.
6. Compare analytic Jacobian blocks with central finite differences.
7. Use mixed absolute/relative gates for near-zero source terms.

### GREEN gate

- Focused C++ tests pass.
- Nonzero analytic/finite-difference Jacobian relative error: `<= 1e-8`.
- Near-zero Jacobian absolute error uses a frozen pre-run threshold.
- Residual and diagnostic replay are identical to `<= 1e-12` relative.
- `git diff --check` and `ascii_sources` pass.
- No unrelated mobility, Poisson, SRH, or carrier-statistics formula changes.

### Commit boundary

Commit tests first where practical, then the minimal opt-in implementation
patch. Do not combine a default change.

## Task 7 - imported-state continuity residual and first QFP update

### Purpose

Address the upstream self-consistent QFP difference that remains after
fixed-state source alignment.

### Actions

At Minimal6 and coarse7x3 exact imported states:

1. Assemble electron and hole continuity residuals using:
   - current production triangle source;
   - opt-in element-edge source;
   - avalanche disabled;
   - SRH disabled control;
   - native-mobility box reconstruction; and
   - unchanged Vela mobility.
2. Record, per node:
   - SG divergence;
   - SRH;
   - electron/hole avalanche source;
   - boundary/Dirichlet contribution;
   - row scaling;
   - normalized and physical residual; and
   - final residual.
3. Compute:
   - carrier-only first Newton update;
   - coupled first Newton update;
   - update direction versus Sentaurus-to-Vela QFP difference; and
   - residual/Jacobian closure.
4. Hold source units, SRH parameters, tolerance, clamp, and continuation fixed.
5. Audit SRH or source-unit scaling only if the term-by-term residual identifies
   it as causal. Symbolic cancellation alone is not evidence for a change.
6. Classify the first bias, node, carrier, and term that drives the QFP branch
   away from Sentaurus.

### Exit gate

- Residual decomposition sums to the assembled residual within `1e-12`
  relative.
- Boundary rows remain unchanged to floating-point precision.
- Jacobian gate from Task 6 remains satisfied.
- First updates exist at `-1`, `-10`, and `-20 V` for both topologies.
- A source-operator change is causal only if it improves the same residual and
  first-update QFP error on both topologies.

### Typed outcomes

- `source_support_causes_qfp_update`;
- `current_coefficient_causes_qfp_update`;
- `boundary_or_contact_model_difference`;
- `srh_or_source_term_difference`;
- `continuation_branch_difference`;
- `proprietary_model_difference`; or
- `operator_improvement_without_qfp_causality`.

## Task 8 - deterministic self-consistent sweeps

### Entry condition

Tasks 6-7 must authorize an opt-in candidate. A fixed-state source improvement
alone is insufficient.

### Actions

1. Rerun the frozen Minimal6 40-state baseline and candidate twice.
2. Run coarse7x3 at exact `-1..-20 V` checkpoints where convergence permits.
3. Preserve every accepted checkpoint and first rejected transition.
4. Compare in dependency order:

   `psi -> QFP -> n/p -> mobility -> directed current -> terminal current -> alpha -> source`.

5. Re-run fixed-state Tasks 3-5 on accepted self-consistent checkpoints.
6. Record:
   - Newton history;
   - residual/source decomposition;
   - terminal-current method closure;
   - KCL;
   - current-growth classification; and
   - deterministic A/B hashes.

### Frozen targets

| Quantity | Target |
|---|---:|
| electrostatic potential maximum | `1e-6 V` |
| electron/hole QFP median/P95 | `0.01/0.025 V` |
| electron/hole density median/P95 | `0.10/0.25 dex` |
| matched-support mobility median/P95 | `0.05/0.20 dex` |
| directed current median/P95 | `0.10/0.25 dex`, 100% nonzero sign |
| total terminal current median | `0.10 dex` |
| active-range impact-source median | `0.30 dex` |
| internal total-current KCL | `1e-8` relative |

### Typed outcomes

- `self_consistent_parity_passed`;
- `operator_improves_but_qfp_misses`;
- `fixed_state_only_improvement`;
- `solver_first_failure`; or
- `model_difference`.

## Task 9 - physical PN2D BV curve and knee comparison

### Entry condition

The coarse self-consistent candidate must not regress QFP, terminal current,
or source conservation. Do not call a Minimal6 curve a physical BV curve.

### Actions

1. Use the fine `reference_tcad/pn2d_sentaurus2018` case.
2. Run baseline and authorized candidate with identical:
   - mesh;
   - contacts;
   - material parameters;
   - bias stepping;
   - continuation;
   - tolerances; and
   - stop conditions.
3. Preserve exact curve rows and first rejected transition.
4. Reuse `scripts/diagnose_pn2d_bv_knee_shape.py`:
   - first one-volt current-growth ratio above `1.5`;
   - first one-volt current-growth ratio above `2.0`; and
   - maximum absolute log-current error over `-20..-10 V`.
5. Reuse `scripts/diagnose_pn2d_bv_knee_gap.py` to distinguish:
   - `latent_turning_point`; and
   - `physics_magnitude_gap`.
6. Add ionization-integral, avalanche source, field maximum, terminal-current,
   and KCL comparisons at each knee marker.
7. Report missing knee markers as `undefined`, never as zero-voltage error.

### Provisional BV gate

- Both Sentaurus knee markers exist in the comparison window.
- Candidate `1.5x` and `2.0x` knee biases are each within `1.0 V` of Sentaurus.
- Maximum log-current error in the knee window: `<= 0.30 dex`.
- No unclassified continuation failure occurs before the Sentaurus knee.
- Terminal-current and KCL method closures remain valid.
- Ionization-integral and source changes have the same causal direction as
  the current knee.

If the reference or candidate has no marker, the outcome is typed
`knee_not_observed`; do not report a numeric voltage difference.

## Task 10 - production decision, validation, review, and commits

### Decision matrix

| Evidence | Allowed action |
|---|---|
| Minimal6 passes but general Tri3 fails | Keep opt-in; no default change. |
| General fixed-state operator passes but first QFP update does not improve | Keep opt-in diagnostic; investigate QFP model. |
| First QFP update improves but self-consistent targets fail | Keep config candidate; no default change. |
| Self-consistent coarse target passes but fine BV knee fails | Keep opt-in; classify physical BV gap. |
| General fixed-state, QFP update, self-consistent, and fine BV gates pass | Propose a separate default-change review. |
| Any improvement requires fitted field, mobility, alpha, or geometry scale | Record model difference; no production patch. |

### Full validation

1. Focused RED/GREEN C++ tests.
2. General-mesh Python regression tests.
3. Two independently generated raw roots.
4. Independent verifier that does not import generator calculations.
5. Minimal6 40-state fixed-state regression.
6. Minimal6 40-state self-consistent regression.
7. coarse7x3 imported-state and self-consistent gates.
8. fine PN2D BV knee audit.
9. Full Release build and CTest.
10. `ascii_sources`, manifest/hash, schema, and scoped-diff checks.

### Independent review

Request:

- scientific review of units, signs, support, contact fallback, source
  identity, QFP causality, and BV-knee claims; and
- code review of configuration propagation, residual/Jacobian consistency,
  general-mesh geometry, tests, and unchanged defaults.

### Commit structure

Keep commits scoped:

1. Task 1 schema and RED tests;
2. Tasks 2-5 remote/export/analyzer/verifier/report;
3. Task 6 opt-in C++ tests and minimal implementation patch, if required;
4. Task 7 imported-state residual and first-update evidence;
5. Task 8 self-consistent candidate configuration and evidence;
6. Task 9 physical BV evidence; and
7. Task 10 decision ledger and review responses.

Never stage unrelated worktree changes or generated simulation roots.

## Formula decision ledger

Maintain one row for each item:

- Poisson electrostatic potential;
- P1 electric-field reconstruction;
- electron/hole QFP gradient;
- Old-Slotboom/BGN carrier statistics;
- low-field mobility;
- high-field mobility;
- SG edge current;
- element-edge current-vector reconstruction;
- Van Overstraeten alpha;
- source mapping and element-vertex measure;
- SRH;
- continuity residual scaling;
- Jacobian;
- nonlinear continuation; and
- BV knee classification.

Each row must be one of:

- `validated_unchanged`;
- `configuration_only`;
- `diagnostic_only`;
- `production_patch`;
- `model_difference`; or
- `insufficient_data`.

No row may inherit another row's conclusion without its own gate.

## Global stop conditions

Stop and preserve partial evidence if:

1. paired branches differ in mesh, TDR, material, model, release, temperature,
   contact bias, or undeclared Math option;
2. topology/cell mapping is not field-verified;
3. a reconstructed edge value is labeled native;
4. contact fallback is cited as the global Sentaurus default;
5. zero support is converted to a finite dex value;
6. low-signal source values dominate an active-range conclusion;
7. carrier sign or edge orientation is ambiguous;
8. source identity or ReadMeasure integration does not close;
9. residual and Jacobian use different mobility, driver, or source mapping;
10. a fitted scale or proprietary inferred value enters a production
    candidate;
11. more than one causal factor changes in an isolated comparison;
12. a production default changes before Task 10;
13. a Minimal6 diagnostic sweep is called a physical BV curve; or
14. a missing BV knee marker is reported as a numeric voltage difference.

## Final success definition

The plan succeeds only when it produces one of two defensible results:

1. a general-mesh, self-consistent, physically validated element-edge
   avalanche operator with an independently reviewed production decision; or
2. a bounded and localized `model_difference` that preserves the accurate
   Vela formulas, identifies the first non-closing QFP/current/source
   dependency, and avoids an unjustified production change.
