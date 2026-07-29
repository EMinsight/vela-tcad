# PN2D avalanche-off SRH spatial-parity follow-up plan

Date: 2026-07-28

Status: proposed for review; execute in a new Codex task after approval.

Starting point: local branch `codex-pn2d-minimal6-operator-audit`, commit
`fa1c343`.

## Background

The corrected two-dimensional continuity-source scale removed the former
two-order-of-magnitude Vela leakage deficit. On the coarse7x3 PN2D reverse
sweep with impact ionization disabled:

| Bias | Vela/Sentaurus current |
|---:|---:|
| -1 V | 3.56 |
| -5 V | 1.14 |
| -10 V | 0.738 |
| -15 V | 0.592 |
| -20 V | 0.835 |

For every nonzero anchor, the independently integrated Vela SRH source agrees
with the summed carrier contact flux to approximately `2e-6` relative or
better. The remaining voltage-shape error is therefore no longer classified
as source loss or false Newton convergence.

The completed `doping_concentration_basis` comparison also establishes:

- `net_doping`, `total_impurity`, and
  `cell_reconstructed_total_impurity` all converge at 21/21 reverse-bias
  points;
- their avalanche-off leakage differences are only at the parts-per-million
  level and have no material physical significance;
- BV keeps `net_doping` as the production baseline;
- forward IV keeps `cell_reconstructed_total_impurity`, where it has a
  separately demonstrated mobility-field and current benefit.

Historical lifetime probes are also frozen evidence: blindly fitting one
constant `taun/taup` value did not explain the complete BV behavior. This plan
must first distinguish state, SRH formula/parameters, spatial integration, and
mesh-resolution effects.

## Primary objective

Identify and correct the dominant cause of the remaining Sentaurus/Vela
avalanche-off BV leakage-shape difference from `-1 V` to `-20 V`, without
regressing continuity closure, forward IV, or avalanche-on infrastructure.

The work must produce one of these evidence-backed classifications:

1. `state_difference`: Vela and Sentaurus carrier/effective-intrinsic-density
   states differ before the SRH operator is applied;
2. `srh_parameter_or_formula_difference`: identical states still produce
   different local SRH rates because the effective lifetime, trap energy, or
   intrinsic-density convention differs;
3. `spatial_support_or_quadrature_difference`: local rates are comparable but
   nodal/cell support, control volumes, or source integration differ;
4. `mesh_resolution_difference`: the discrepancy decreases systematically
   under junction-focused mesh refinement;
5. `mixed`: more than one independently quantified contribution is required.

## Frozen production configuration

During Tasks 1-4, hold these choices fixed:

- Vela BV mobility basis: `net_doping`;
- Vela IV mobility basis: `cell_reconstructed_total_impurity`;
- impact ionization: disabled;
- SRH and Old Slotboom BGN: enabled exactly as in the aligned BV decks;
- temperature, contacts, device depth, current normalization, solver
  tolerances, bias points, and initialization: unchanged;
- sweep points: `0, -1, ..., -20 V`;
- coarse7x3 mesh and doping profile: unchanged unless a task explicitly
  belongs to the mesh-refinement matrix.

Do not alter avalanche coefficients, high-field mobility, contact models, or
global production defaults in this plan.

## Common metrics and artifact contract

Every experiment must write a machine-readable manifest containing:

- Git commit and dirty-state flag;
- Vela and Sentaurus input/config hashes;
- mesh, doping, material, and state hashes;
- bias, model switches, lifetime/trap parameters, and unit system;
- convergence status and continuity-closure status;
- terminal electron, hole, and total currents;
- integrated electron and hole SRH source;
- wall-clock time.

At each bias, report:

- signed and absolute terminal current;
- `log10(|I_Vela/I_Sentaurus|)`;
- electron and hole global closure error;
- total-current terminal closure;
- integrated positive generation and negative recombination separately;
- SRH source centroid and the 10%, 50%, and 90% cumulative-source positions
  along the junction-normal coordinate;
- depletion-support width under one fixed, documented depletion indicator.

Use `1e-30 A/um` only as a plotting/logarithm floor. It must not modify the
simulation or linear/log error calculations above the floor.

## Task 1 - Freeze the baseline and add spatial SRH diagnostics

### Goal

Create one reproducible baseline dataset that shows where Vela and Sentaurus
generate avalanche-off leakage, not only how much terminal current they
produce.

### Implementation steps

1. Add a runner derived from the existing avalanche-off comparison scripts.
   It must generate the aligned 21-point Vela and Sentaurus datasets without
   changing either deck.
2. Extend Vela diagnostics to export, per node and per triangle:
   - coordinates and control-volume/cell area;
   - `psi`, `n`, `p`, effective intrinsic density, and net doping;
   - raw SRH rate;
   - signed and absolute area-integrated SRH contribution;
   - depletion indicator and distance from the metallurgical junction.
3. Import the equivalent available Sentaurus fields from TDR/CSV. Mark every
   unavailable native field explicitly; do not silently reconstruct it.
4. Build cumulative-source profiles along the junction-normal coordinate and
   generate:
   - terminal avalanche-off IV comparison;
   - local SRH-rate comparison;
   - cumulative SRH-source comparison;
   - source centroid and percentile-width table.
5. Add a focused test for area-weighted source integration using a synthetic
   triangular mesh with a known constant and linear source field.

### Deliverables

- `scripts/run_pn2d_bv_off_srh_spatial_audit.py`
- `scripts/build_pn2d_bv_off_srh_spatial_report.py`
- a dated report under `docs/validation/`
- CSV/JSON data and figures under an ignored `build-release/` report directory

### Acceptance criteria

- Vela and Sentaurus both converge at 21/21 points.
- Vela electron and hole source/contact closure is `<= 1e-5` relative at every
  nonzero bias.
- Vela total-current terminal closure is `<= 1e-20 A/um` absolute.
- Reintegrating exported Vela nodal/cell SRH data reproduces the solver's
  reported integrated SRH source within `1e-6` relative.
- The synthetic constant/linear source tests pass to `1e-12` relative.
- Every plotted value is traceable to a manifest and CSV column.

### Stop condition

If exported Vela source does not reproduce the solver integral, stop before
cross-simulator interpretation and repair the diagnostic/integration path.

## Task 2 - Perform a same-state SRH operator decomposition

### Goal

Separate state differences from SRH formula and spatial-integration
differences.

### Implementation steps

1. At `-1, -5, -10, -15, -20 V`, import the exact Sentaurus
   `psi/n/p` fields and all available intrinsic-density/BGN fields onto the
   identical mesh.
2. Evaluate the Vela SRH operator on:
   - the self-consistent Vela state;
   - the imported Sentaurus state;
   - a hybrid state changing only `n,p`;
   - a hybrid state changing only effective intrinsic density/BGN inputs.
3. For every state, integrate the same local rate using:
   - Vela nodal control volumes;
   - element-based triangle quadrature;
   - any recoverable Sentaurus native source integral/support.
4. Decompose the log-current gap into:
   - state contribution;
   - effective-intrinsic-density/BGN contribution;
   - local SRH formula/parameter contribution;
   - integration/support contribution.
5. Add conservation and permutation tests for the element integration path.

### Acceptance criteria

- Imported fields cover 100% of required semiconductor nodes, with no duplicate
  node IDs and no implicit zero fill.
- Same-state repeated runs reproduce integrated source within `1e-10`
  relative.
- Nodal and element integration agree within `1%` on the coarse mesh or the
  report identifies and spatially localizes the cells responsible for a
  larger difference.
- At least 90% of the terminal log-current gap at each anchor is assigned to
  named terms, or the residual unassigned term is reported as a blocker.
- No production formula is changed during this diagnostic task.

### Decision gate

- If imported Sentaurus state plus the existing Vela SRH operator matches the
  Sentaurus integral within `0.05 dex`, classify the primary issue as
  `state_difference`.
- If local rates disagree before integration by more than `0.05 dex`, proceed
  to Task 3.
- If local rates agree but integrated sources differ by more than `0.05 dex`,
  classify the issue as `spatial_support_or_quadrature_difference` and carry
  that evidence into Task 4.

## Task 3 - Audit Sentaurus SRH parameters and effective lifetime

### Goal

Determine whether Sentaurus uses a materially different SRH formula or
effective parameter set, without treating free lifetime fitting as proof.

### Implementation steps

1. Seal the exact Sentaurus BV command file, material parameter source,
   release, log, and all active Physics/Math settings.
2. Build a parameter-parity table for:
   - `taun` and `taup`;
   - trap energy/reference level;
   - electron and hole capture conventions;
   - effective intrinsic density and BGN coupling;
   - any field-, doping-, or temperature-dependent lifetime switches.
3. Where Sentaurus exports enough fields, invert the SRH equation to infer an
   effective lifetime or denominator at each supported node. Report its
   variation with bias, doping, electric field, carrier density, and position.
4. Run a preregistered, minimal candidate matrix only:
   - exact current Vela parameters;
   - exact documented/deck-resolved Sentaurus parameters;
   - at most one evidence-derived candidate if the inversion rejects both.
5. Use `-1, -10, -20 V` as calibration/diagnostic anchors and
   `-5, -15 V` as hold-out validation points.

### Acceptance criteria

- Every candidate parameter has a deck, manual, material-file, log, or
  inversion provenance; no unexplained fitted constant is accepted.
- A constant-lifetime claim is allowed only if the inferred effective value
  has `<= 5%` robust spread over the source-dominant region at all anchors.
- An evidence-derived candidate must improve hold-out log-current RMSE by at
  least 30% relative to the frozen baseline.
- It must not worsen any calibration or hold-out anchor by more than
  `0.05 dex`.
- Continuity closure remains `<= 1e-5` relative.

### Stop condition

If no provenance-backed candidate passes the hold-out criteria, keep the
existing SRH parameters and classify parameter tuning as rejected.

## Task 4 - Run a junction-focused mesh-convergence study

### Goal

Measure whether coarse spatial support is responsible for the leakage-shape
difference and determine the resolution needed for a mesh-independent answer.

### Implementation steps

1. Define at least three nested mesh levels:
   - `M0`: existing coarse7x3;
   - `M1`: approximately 2x finer spacing in the depletion/junction region;
   - `M2`: approximately 4x finer spacing in that region.
2. Preserve geometry, doping profile, contacts, material interfaces, device
   depth, and physical settings. Record mesh-quality statistics and hashes.
3. Use the same physical mesh in Vela and Sentaurus at each level. A
   simulator-specific remesh is not a valid paired comparison.
4. Run avalanche-off at all 21 bias points on each mesh.
5. Compare current, integrated source, centroid, percentile width, peak rate,
   and depletion width across mesh levels.
6. Estimate observed mesh-convergence order where three valid levels exist.

### Acceptance criteria

- Both simulators converge at 21/21 points on every accepted mesh.
- Doping dose per device depth changes by `< 0.1%` between mesh levels.
- Contact locations and device dimensions are identical within the mesh input
  precision.
- On the two finest levels, integrated SRH source changes by `< 2%` at each
  anchor for a mesh-independent result.
- A mesh-root-cause claim requires the Vela/Sentaurus log-current RMSE to
  improve monotonically and by at least 30% from M0 to M2.
- Closure remains `<= 1e-5` relative at every nonzero bias.

### Decision gate

- If both simulators converge to the same curve within `0.05 dex`, adopt a
  documented PN2D junction-resolution requirement rather than changing SRH
  physics.
- If each simulator is mesh-converged but their difference remains above
  `0.05 dex`, return to the Task 2 decomposition and implement only the
  identified operator/parameter difference.

## Task 5 - Implement the minimal evidence-selected correction

### Goal

Apply one narrowly scoped correction supported by Tasks 2-4.

### Allowed correction classes

- SRH parameter/configuration support with explicit JSON schema and defaults;
- effective-intrinsic-density/BGN coupling correction;
- nodal/element SRH source-support or quadrature correction;
- mesh-generation/refinement requirement and validation tooling;
- a combination only when the decomposition quantitatively proves multiple
  independent causes.

### Implementation steps

1. Write a failing focused C++ or Python regression reproducing the selected
   discrepancy.
2. Implement the smallest production change that makes the focused regression
   pass.
3. Keep backward-compatible defaults unless the parity evidence explicitly
   proves the old default physically incorrect.
4. Update configuration schema, templates, and diagnostics as applicable.
5. Re-run Tasks 1-4 metrics with the candidate enabled.

### Acceptance criteria

- Focused RED test fails before and passes after the correction.
- Avalanche-off median absolute log-current error over 21 points is
  `<= 0.05 dex`; maximum error is `<= 0.10 dex`.
- If that absolute target is not reached, the correction must still reduce
  21-point log-current RMSE by at least 50% and leave a quantified residual
  follow-up.
- Electron and hole closure remains `<= 1e-5` relative.
- All existing SRH, SG flux, Newton, DC sweep, scaling, and template tests pass.

## Task 6 - Production qualification and final decision

### Goal

Verify that the BV fix is physical, numerically stable, and isolated from the
established IV configuration.

### Validation matrix

1. BV avalanche-off, `0..-20 V`, 21 points, `net_doping`.
2. BV avalanche-on production baseline, `0..-20 V`, 21 points,
   `net_doping`.
3. Forward IV, `0..20 V`, 201 points,
   `cell_reconstructed_total_impurity`.
4. Focused same-state SRH audit at `-1, -5, -10, -15, -20 V`.
5. Full Debug CTest suite.

### Acceptance criteria

- All required bias points converge; no point may be silently interpolated.
- Avalanche-off meets the Task 5 current-error and closure gates.
- Forward-IV current changes by `< 0.5%` at every nonzero reference anchor and
  retains 201/201 convergence.
- Avalanche-on must not move by more than `0.05 dex` below its established
  multiplication-onset region unless the change is directly predicted by the
  corrected SRH baseline.
- Full CTest passes with zero failures.
- The final report states separately:
  - whether the scientific root cause is closed;
  - whether a production-default change is authorized;
  - the final BV and IV `doping_concentration_basis` choices.

## Execution order and commit boundaries

Execute strictly in this order:

1. Task 1 baseline and diagnostics;
2. Task 2 same-state decomposition;
3. Task 3 parameter audit only if its decision gate is entered;
4. Task 4 mesh convergence;
5. Task 5 evidence-selected fix;
6. Task 6 qualification.

Recommended local commits:

1. `test: add PN2D SRH spatial diagnostics`
2. `test: decompose PN2D same-state SRH source`
3. `test: add PN2D SRH parameter and mesh matrix`
4. `fix: align PN2D avalanche-off SRH leakage`
5. `docs: record PN2D SRH spatial-parity decision`

Do not commit generated simulation outputs. Do not stage `tmp/`.

## New-task handoff prompt

Use the following prompt in the new Codex task:

```text
请执行
docs/superpowers/plans/2026-07-28-pn2d-bv-off-srh-spatial-parity.md。
从任务1开始，严格遵守各任务的验收条件和停止条件。先确认当前分支包含
fa1c343，并保持 BV 使用 net_doping、正向 IV 使用
cell_reconstructed_total_impurity。每完成一个决策门请汇报证据，再继续下一任务；
不要提交 generated simulation outputs，也不要暂存 tmp/。
```
