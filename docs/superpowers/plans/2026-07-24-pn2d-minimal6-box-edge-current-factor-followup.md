# PN2D Minimal6 box-edge current-factor follow-up plan

Date: 2026-07-24

Status: Tasks 1-10 complete; validation, review requests, and scoped commits
finished.

## Execution summary

- Tasks 1-4 sealed the native-element electric-field replay and box-current
  projection. The mobility/current replay passes magnitude and sign gates but
  the fixed-state KCL gate remains a typed bounded failure.
- Task 5 classified the three-factor result as `interaction_dominant`.
- Task 6 established causal first-update improvement for the existing global
  `electric_field` mobility configuration.
- Task 7 authorized a config-only candidate, not a formula patch.
- Task 8 classified the 40-state candidate as
  `config_improves_but_misses_target`; production defaults remain unchanged.
- Task 9 corrected two diagnostic source contracts and classified the
  production triangle impact-source difference as
  `current_support_dominant`.
- No fitted mobility, field, saturation, impact, or geometry coefficient was
  introduced.


## Relationship to the existing plan

This plan supersedes only the `Immediate next action` in
`2026-07-24-pn2d-minimal6-current-alignment-phase-b-g-revised.md`.
Completed Phase B-F evidence remains frozen. Final Phase G production
acceptance is deferred until the box-edge mobility/current experiment in this
plan is complete.

## Objective

Determine, without parameter fitting, how much of the remaining fixed-state
and self-consistent current difference is caused by:

1. the low-field mobility coefficient;
2. the high-field driving-force definition;
3. native element versus production global-edge support;
4. the resulting continuity residual and first nonlinear update; and
5. the downstream self-consistent QFP branch.

The immediate target is a diagnostic-only box-edge branch:

`Sentaurus native low-field element mobility`
`+ native element electric field`
`+ documented high-field law`
`+ coefficient-weighted box-edge projection`.

No production mobility, SG, QFP, SRH, impact, or Poisson formula may be
changed before the fixed-state and first-update gates below pass.

## Frozen evidence

The following observations are inputs, not values to be refitted:

- exact lattice: `mirror/sketch x -1..-20 V`;
- 40 exact states, 160 elements, 320 carrier-element samples, and 400 active
  carrier-edge samples;
- imported Sentaurus potential/QFP plus Vela BGN recomputes carrier density
  within `4.426181e-6 dex`;
- imported-state Vela-mobility current median error:
  `0.062727 dex` electron and `0.057626 dex` hole;
- self-consistent directed-current median error:
  `0.823474 dex` electron and `0.939922 dex` hole;
- paired imported-state improvement:
  `0.721335 dex` electron and `0.848517 dex` hole;
- all 80 self-consistent sign mismatches are on central edge `1-5`; the
  corresponding reference currents are approximately `1e-20 A/um`;
- box geometry contributes zero error;
- the native Sentaurus-mobility box reconstruction closes total terminal
  current to `1.137692e-7` relative and internal total-current KCL to
  `3.163255e-9`;
- Sentaurus native low-field versus Vela cell-average low-field mobility
  median errors are `0.139859 dex` electron and `0.085160 dex` hole;
- using Sentaurus native low-field mobility and native element electric field
  replays final element mobility over all 40 states with:
  - electron median/P95/maximum
    `0.001028/0.019657/0.036146 dex`;
  - hole median/P95/maximum
    `0.000496/0.017325/0.027579 dex`;
- the corresponding exported-native-QFP-gradient replay errors are
  `0.497771/0.618552/0.641861 dex` for electrons and
  `0.023278/0.071141/0.098285 dex` for holes; and
- Sentaurus directed-edge current remains unavailable. Every Sentaurus edge
  value must remain labeled `box_operator_reconstruction`.

## Causal model

The experiment must preserve the following dependency order:

`psi and QFP state`
`-> carrier density`
`-> low-field mobility`
`-> high-field drive`
`-> element final mobility`
`-> box-edge mobility`
`-> SG edge current`
`-> node continuity residual`
`-> first nonlinear update`
`-> self-consistent QFP state`
`-> terminal current and impact source`.

A downstream improvement is not evidence for an upstream formula change
unless every earlier dependency is held fixed or explicitly replaced.

## Work products

Expected new diagnostic files:

- `scripts/pn2d_minimal6_diagnostics/highfield_box_replay.py`;
- `scripts/diagnose_pn2d_minimal6_highfield_box_current.py`;
- `scripts/verify_pn2d_minimal6_highfield_box_current.py`;
- `tests/regression/test_pn2d_minimal6_highfield_box_current.py`;
- `docs/validation/pn2d_minimal6_highfield_box_current_2026-07-24.md`; and
- two independently generated evidence roots under `build-release/`.

Production files are out of scope until Task 7 authorizes a candidate.

## Task 1 - seal the experiment contract and create the RED test

### Purpose

Prevent a later branch from silently mixing node, element, global-edge, or
topology-cell order.

### Actions

1. Add a regression test that requests the new typed branch
   `sentaurus_lowfield_element_electric_field`.
2. Require the diagnostic schema to record:
   - topology, bias, carrier, Vela triangle id, Sentaurus region-cell id;
   - low-field source and SHA-256;
   - electric-field source and SHA-256;
   - high-field exponent and saturation velocity;
   - element-edge box coefficient;
   - zero/reference-missing status; and
   - reconstruction label.
3. Require the mapping:
   - mirror: `0->0, 1->1, 2->2, 3->3`;
   - sketch: `0->0, 1->3, 2->2, 3->1`.
4. Make the initial test fail because the branch and output schema do not yet
   exist.

### RED gate

- The focused regression test fails for the missing typed branch.
- It must not fail due to missing historical evidence, locale, newline, or
  platform-dependent path formatting.

## Task 2 - validate native low-field mobility reuse

### Purpose

Verify the assumption that one HighFieldSaturation-off result per topology can
be reused over all biases.

### Actions

1. Retain the completed `-20 V` mirror and sketch controls.
2. Regenerate HighFieldSaturation-off controls at `-1 V` and `-10 V` for both
   topologies.
3. Keep `DopingDependence`, temperature, mesh, material parameters, SRH,
   avalanche, Old Slotboom, and contact definitions unchanged.
4. Export native element:
   - electron/hole mobility;
   - electric field;
   - electron/hole QFP gradient; and
   - current-density vectors as a non-authoritative control.
5. Compare the four element low-field mobility values at `-1`, `-10`, and
   `-20 V`.

### Exit gate

- Six exact control runs pass.
- Maximum relative low-field mobility variation across bias is at most
  `1e-10`.
- Mirror/sketch differences are explained exactly by the verified cell
  permutation.
- If this gate fails, stop reuse and generate the full 40-state low-field
  lattice before continuing.

## Task 3 - GREEN native-element high-field replay

### Purpose

Turn the existing exploratory inversion into a sealed, independently verified
element-level reconstruction.

### Actions

1. Compute for every carrier element:

   `mu_E = mu0_sent / (1 + (mu0_sent * |E_element| / vsat)^beta)^(1/beta)`.

2. Retain comparison branches:
   - exported native QFP gradient;
   - affine triangle QFP gradient;
   - native element electric field;
   - low-field-only; and
   - native exported final mobility reference.
3. Report mobility error by carrier, topology, bias, cell, and current
   weighting.
4. Preserve inversion results separately from direct final-mobility replay.
5. Run the independent verifier and generate byte-identical A/B roots.

### Exit gate

For the electric-field branch:

- median mobility error at most `0.005 dex`;
- P95 at most `0.03 dex`;
- maximum at most `0.05 dex`;
- it beats the triangle-QFP branch on both median and P95 for both carriers;
- 320 samples are valid and finite; and
- A/B roots are byte-identical.

Failure produces a typed `high_field_replay_mismatch`; it does not authorize
parameter fitting.

## Task 4 - project the candidate to box edges and replay current

### Purpose

Measure whether the element-level mobility improvement survives the exact
box-edge projection and SG current operator.

### Actions

1. For each active global edge compute:

   `mu_box = sum(kappa_cell_edge * mu_cell) / sum(kappa_cell_edge * active)`.

2. Build the following edge-mobility branches:
   - native Sentaurus final element mobility reference;
   - Sentaurus low-field plus native electric field;
   - Sentaurus low-field plus triangle QFP gradient;
   - Sentaurus low-field plus exported native QFP gradient;
   - Vela low-field plus native electric field;
   - unchanged Vela production mobility; and
   - constant-mobility control.
3. Hold `psi`, QFP, density, geometry, carrier sign, edge orientation, and SG
   formula identical across all branches.
4. Replay all 40 states and record:
   - coefficient-weighted edge mobility;
   - directed electron/hole current;
   - sign agreement;
   - absolute and current-weighted dex error;
   - contact carrier current;
   - total terminal current;
   - internal total-current KCL; and
   - geometric-zero edge status.
5. Compare every candidate first with the native-mobility reconstructed edge
   reference, then with Sentaurus terminal current.

### Exit gate

For the low-field-plus-electric-field branch relative to the native-mobility
box reference:

- 400 active carrier edges are present;
- median directed-current error at most `0.01 dex`;
- P95 at most `0.05 dex`;
- sign agreement is 100%;
- zero-coefficient edges remain exact typed zeros;
- total-terminal difference is at most `2%`;
- internal total-current KCL remains below `1e-8`; and
- the original native-mobility reference still closes Sentaurus terminal
  current below `2e-7`.

Report the central `1-5` edge both in dex and absolute `A/um`. It must not
dominate a current-weighted conclusion merely because its reference is near
zero.

## Task 5 - perform a factorial current attribution

### Purpose

Separate low-field coefficient, high-field drive, and support effects without
making the conclusion depend on replacement order.

### Actions

1. Evaluate the fixed imported state over the factorial controls:
   - low-field source: Vela or Sentaurus;
   - drive: global-edge QFP, triangle QFP, or native element electric field;
   - support: one global edge or coefficient-weighted native elements.
2. Compute all admissible replacement orders.
3. Report:
   - paired incremental contribution for every order;
   - order-averaged Shapley contribution;
   - interaction term;
   - unweighted median/P95/max; and
   - absolute-reference-current-weighted contribution.
4. Separate ordinary active edges from the central `1-5` tail.
5. Require signed-log decomposition closure for every sample.

### Exit gate

- Maximum decomposition closure error is at most `1e-12 dex`.
- Every reported contribution has a same-edge paired baseline.
- The ranking is stable for unweighted and current-weighted results, or the
  disagreement is explicitly classified.
- No fitted mobility, field scale, or edge coefficient is introduced.

Typed outcomes:

- `low_field_coefficient_dominant`;
- `driving_force_dominant`;
- `support_dominant`;
- `interaction_dominant`; or
- `mixed_bounded_contributions`.

## Task 6 - carry the branch into the continuity residual

### Purpose

Determine whether the improved current coefficient changes the imported-state
continuity residual and first QFP update enough to explain the remaining
self-consistent branch difference.

### Actions

1. At all 40 imported Sentaurus states, assemble node residuals with:
   - unchanged Vela production mobility;
   - native Sentaurus final element mobility;
   - Sentaurus-low-field plus electric-field replay;
   - triangle-QFP replay; and
   - constant mobility.
2. Preserve separate terms:
   - SG divergence;
   - SRH;
   - impact generation;
   - boundary/Dirichlet row;
   - normalization; and
   - final residual.
3. At `-1`, `-10`, and `-20 V`, compute the first carrier-only and coupled
   nonlinear updates for each branch.
4. Do not combine a mobility change with a source-unit, tolerance,
   continuation, or clamp change.
5. Re-run analytic versus finite-difference Jacobian checks if a branch is
   moved into C++ diagnostics.

### Exit gate

- Edge-to-node residual replay closes the recorded operator below `1e-12`
  relative where the reference is nonzero.
- Boundary rows remain unchanged to floating-point roundoff.
- Analytic/finite-difference Jacobian block difference remains below `1e-8`.
- The first-update comparison is available for both carriers and both
  topologies at all three selected biases.
- A mobility candidate is causal only if it reduces the same imported-state
  residual and first-update QFP error on both topologies without a new sign or
  convergence regression.

If element current improves but the first update does not, classify the result
as `current_coefficient_improvement_without_qfp_causality`.

## Task 7 - production-decision gate

### Purpose

Choose the smallest justified action only after Tasks 1-6.

### Decision matrix

| Evidence | Allowed action |
|---|---|
| Electric-field replay fails element or edge gates | Keep diagnostic evidence; no mobility change. |
| Replay passes, and an existing `electric_field` configuration reproduces it on the required support | Change only the comparison configuration; no formula patch. |
| Replay passes, but the difference is solely element-to-global-edge support | Propose a discretization design separately; do not hide it in a field-scale change. |
| Fixed-state residual and first update improve causally with an existing configuration | Authorize a config-only self-consistent candidate sweep. |
| A Vela implementation defect is independently reproduced in residual and derivative | Add a minimal RED/GREEN production patch limited to that defect. |
| Improvement requires fitted scales or undocumented vendor behavior | Record `model_difference`; no production patch. |

### Required ledger

Create a formula-decision ledger for:

- carrier-density transformation;
- low-field mobility;
- high-field saturation law;
- high-field driving force;
- element-to-edge projection;
- SG current;
- SRH source;
- impact source;
- Poisson; and
- nonlinear continuation.

Each row must state `unchanged`, `configuration-only`, `diagnostic-only`, or
`production-patch`, with an evidence link.

## Task 8 - config-only self-consistent candidate sweep

### Entry condition

Task 7 must authorize an existing configuration, preferably
`high_field_driving_force: electric_field`. Do not implement a new production
operator in this task.

### Actions

1. Run baseline and candidate Release sweeps on
   `mirror/sketch x -1..-20 V`.
2. Preserve exact checkpoints, first rejected transition, solver histories,
   and branch classification.
3. Compare in dependency order:

   `psi -> QFP -> density -> mobility -> edge current -> terminal current -> impact`.

4. Generate two independent roots per branch.
5. Re-run the Phase B/C fixed-state gates unchanged.

### Frozen targets

| Quantity | Target |
|---|---:|
| electrostatic potential maximum | `1e-6 V` |
| electron/hole QFP median/P95 | `0.01/0.025 V` |
| electron/hole density median/P95 | `0.10/0.25 dex` |
| directed current median/P95 | `0.10/0.25 dex`, 100% sign |
| total terminal current median | `0.10 dex` |
| impact-source median | `0.30 dex` |

The candidate must also improve both carrier QFP metrics relative to the
current Phase F baseline. A downstream current improvement with worse QFP is
not sufficient.

Typed outcomes:

- `config_parity_passed`;
- `config_improves_but_misses_target`;
- `config_no_causal_improvement`;
- `solver_first_failure`; or
- `model_difference`.

## Task 9 - impact and avalanche follow-up

### Entry condition

Do not start until QFP and current alignment are classified by Task 8.

### Actions

1. Recompute electric field, QFP gradient, carrier current, alpha, and
   integrated avalanche source on identical state/support.
2. Separate mobility-drive changes from impact-ionization-drive changes.
3. Retain the zero-area diagonal-edge geometry contract.
4. Compare source-to-node mapping and Jacobian independently.

### Exit gate

- Current support and sign close before alpha is evaluated.
- Alpha-law inputs are identical before formula parity is claimed.
- Source integration and node mapping close independently.
- No impact coefficient is fitted to compensate for a QFP/current error.

## Task 10 - independent validation, review, and commit

### Validation

1. Focused RED/GREEN regression tests.
2. A/B byte-identical evidence roots.
3. Independent verifier on raw samples and summaries.
4. Phase B terminal/KCL replay.
5. Phase C 40-state staged replacement replay.
6. Phase E residual/Jacobian guards.
7. Phase F baseline comparison.
8. Full Release build and CTest.
9. `ascii_sources`, schema, hash, and scoped-diff checks.

### Review

Request two independent reviews:

- scientific review of units, signs, support, causal claims, and stop
  conditions; and
- code review of every changed production or configuration line.

### Commit structure

Keep commits scoped:

1. diagnostic schema, RED/GREEN tests, replay, verifier, and evidence report;
2. config-only candidate, if authorized;
3. production patch, only if separately authorized by Task 7; and
4. final validation ledger and review responses.

Do not stage unrelated existing worktree changes.

## Independent source-unit guard

The mobility experiment must hold the current source-unit implementation
fixed. The user's concern that an SRH scale may cancel in a normalized
continuity equation remains a separate audit question and must not be
answered by mobility-current improvement.

Before final production acceptance:

1. compare the historical forward-IV regression with and without the
   source-unit change;
2. write the dimensional equation before and after row normalization;
3. demonstrate whether the scale cancels in residual, Jacobian, update, and
   terminal current, rather than in the symbolic source term alone; and
4. retain or revert that patch based on its own RED/GREEN evidence.

No mobility result may be used to justify an SRH/source-unit decision, and no
source-unit result may be used to claim mobility-drive parity.

## Global stop conditions

Stop and preserve partial evidence if:

1. topology, bias, source hash, cell mapping, units, or carrier orientation
   differs across paired branches;
2. low-field mobility is not bias invariant and the full low-field lattice is
   unavailable;
3. Vela production current cannot be replayed below `1e-12` relative;
4. the native-mobility Sentaurus reconstruction loses its terminal/KCL gates;
5. a geometric zero is converted to a finite dex value;
6. a native QFP-gradient field is relabeled as an internal Sentaurus drive;
7. a fitted scale, saturation velocity, beta, or edge coefficient is
   introduced;
8. a production change is made before fixed-state residual and first-update
   causality pass;
9. mirror/sketch disagree after the verified region-cell permutation without
   a localized reason; or
10. a candidate changes more than one causal factor in the same comparison.

## Final execution result

Tasks 1-10 are complete. The `electric_field` high-field mobility candidate
improves the QFP and current comparisons but misses the frozen acceptance
targets and slightly worsens the production triangle impact-source comparison.
It therefore remains an explicit diagnostic candidate rather than the
production default.

The next experiment is outside this completed plan: align the production
triangle current proxy with the independently closed box-current support before
considering any impact-ionization formula or default-configuration change.
