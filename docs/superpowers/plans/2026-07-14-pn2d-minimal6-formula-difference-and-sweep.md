# PN2D Minimal6 Formula-Difference And Diagnostic-Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce reproducible fixed-state causal evidence for the Sentaurus/Vela avalanche-source gap, then test whether the same factors remain dominant in independent sketch/mirror minimal6 nonlinear diagnostic sweeps without changing production physics or solver defaults.

**Architecture:** First recover or regenerate authoritative Sentaurus inputs and fail closed on provenance. Put units, topology, formula reconstruction, counterfactual evaluation, schemas, and plotting behind a small reusable Python package; keep the fixed-state and sweep CLIs as orchestration layers. Treat fixed-state diagnosis and nonlinear sweep as separate publishable phases joined only through versioned ledger records and explicit state identities.

**Tech Stack:** Python 3 standard library plus existing NumPy/Matplotlib/Pillow dependencies, C++20/CMake/Ninja, Catch2, MSYS2 UCRT64, Sentaurus Device 2018 through the existing `sentaurus` SSH workflow.

## Technical Review And Revisions

1. The historical handoff names `minimal6_states_live_20260713_v2`, but that ignored directory is not currently reproducible from the tracked tree. Historical row counts and errors are context, not a current data gate. Phase 0 therefore precedes all report work.
2. `AvalancheGeneration` must be the external source anchor. Sentaurus `alpha * |J| / q` is a separately labelled reconstruction and cannot close the native-source residual by definition.
3. A single ordered waterfall is path-dependent. The report must emit forward and reverse paths, and emit a 2x2 matrix for adjacent factors whose path contributions differ by more than `0.3 dex`.
4. BGN, mobility, current semantics, driving field, alpha law, geometric volume, and node mapping do not all live on the same support. Every substitution must pass through an explicitly named support-conversion operator; direct node/edge/triangle mixing is forbidden.
5. Sentaurus ionization integrals are path-dependent diagnostics, not local alpha or local generation. Store them in the ledger, but do not use them as interchangeable waterfall factors.
6. Do not grow `scripts/audit_pn2d_minimal6_fixed_state.py` into a second analysis framework. Preserve it as the imported-state parity/replay gate and build the new causal analysis from focused modules.
7. Vela's current `node_doping_file` accepts net doping. OldSlotboom comparison needs raw donor and acceptor concentrations as well as net doping; the sweep manifest must record which representation each solver consumed.
8. Integer checkpoint orchestration must not imply that Vela internally used 1 V steps. The driver may use adaptive continuation inside each `[k, k-1] V` segment, but only a solver-accepted endpoint exactly equal to the integer target may be published.
9. The two new public schemas describe reports, not raw solver states. Reuse `vela.pn2d_minimal6_states.v1` for authoritative fixed-state inputs and give Sentaurus sweep states a separate internal manifest identity.

## Global Constraints

- Do not modify production C++ APIs, production formulas, default model parameters, source mapping, or solver defaults.
- Do not commit generated Sentaurus, Vela, plot, or report artifacts; keep them below ignored `build-release/` roots.
- Use `D:\msys64\ucrt64\bin` and `D:\msys64\usr\bin` first on `PATH`.
- Required fixed-state matrix is exactly `sketch/mirror x 0/-12/-19 V`, with `6 nodes / 9 unique edges / 4 CCW triangles` per state.
- Required fixed-state row counts remain exactly `36 node / 54 edge / 24 triangle`.
- Imported-state parity must be `<1e-12`; same-formula C++/Python hybrid error must be `<5e-12`.
- Bias identity tolerance is `1e-12 V`; nearest-state and interpolated-state substitution are forbidden.
- Geometry support `<=1e-27 m^2` or both compared source magnitudes `<=1e-285` is classified as `geometric_zero`; it must not receive an artificial log ratio.
- `0 V` is a zero-field/numerical-floor control and is excluded from dominant-root-cause scoring.
- A factor is dominant only if its signed direction agrees at `-12 V` and `-19 V` for both topologies and it explains at least 50% of the absolute native-source log gap under the report's stated path convention.
- Every external log gap must equal named factor contributions plus a named `sentaurus_internal_semantics_residual`, within the report serialization tolerance.
- Every curve and report must say `minimal6 diagnostic sweep; not a physical BV curve`.

---

### Task 1: Authoritative-State Recovery And Extended Field Contract

**Files:**
- Modify: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd`
- Modify: `scripts/export_pn2d_minimal6_states.py`
- Modify: `tests/regression/test_pn2d_minimal6_state_export.py`
- Create: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`

**Interfaces:**
- Consumes: an archived `minimal6_states_live_20260713_v2` directory or the existing controlled Sentaurus export path.
- Produces: a validated `vela.pn2d_minimal6_states.v1` root with the six exact states and extended required-field entries.

- [ ] **Step 1: Write RED contract tests**

Add cases requiring scalar `AvalancheGeneration`, vector `eVelocity`/`hVelocity`, and scalar electron/hole/mean ionization-integral exports with explicit units and complete global-node mapping. Add rejection cases for missing fields, wrong units/components, an inexact requested/actual bias, and an archive whose manifest or member hash changed.

- [ ] **Step 2: Extend the SDevice plot contract**

Add the exact Sentaurus plot quantities supported by the installed 2018 deck syntax for native avalanche generation, carrier velocity, and carrier/mean ionization integrals. Keep the existing potential, QF, density, electric field, current, mobility, and alpha fields unchanged. If the installed version exposes an integral only under a different canonical name, record the raw name and normalized ledger name in `field_manifest.json`; do not synthesize the field.

- [ ] **Step 3: Extend `REQUIRED_FIELDS` and validation**

Represent every required field as `(normalized_name, raw_name, components, unit, semantic_role)`. Require exactly one matching region-0 field and reject duplicate aliases. Preserve raw field filenames and hashes in each state record.

- [ ] **Step 4: Attempt archive recovery before remote execution**

Search only declared local archive/checkpoint locations for the exact run ID. Validate the original manifest hash, all member hashes, six-state identity, exact biases, topology, fields, and `outputs_complete=true`. Never rewrite the recovered manifest.

- [ ] **Step 5: Regenerate under a new run ID if recovery fails**

Run `scripts/export_pn2d_minimal6_states.py` through the existing `sentaurus` SSH target with a new timestamped run ID. A network/authentication failure stops the plan here after preserving the partial manifest and logs.

- [ ] **Step 6: Replay the C++ audit and seal provenance**

For every state, generate topology-matched `mesh.json` and immutable `audit.json`, run `build-release/pn2d_minimal6_operator_audit.exe`, and record executable/config/input/output SHA-256 values and the full argv array. Then run the existing fixed-state CLI and require `36/54/24`, parity `<1e-12`, and formula error `<5e-12`.

- [ ] **Step 7: Verify focused regressions**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_fixed_state_audit -v
```

Expected: all tests pass; the validation document names either the recovered immutable run ID or the newly generated run ID.

### Task 2: Shared Quantity-Ledger Core And Report Schemas

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/__init__.py`
- Create: `scripts/pn2d_minimal6_diagnostics/contracts.py`
- Create: `scripts/pn2d_minimal6_diagnostics/units.py`
- Create: `scripts/pn2d_minimal6_diagnostics/ledger.py`
- Create: `scripts/pn2d_minimal6_diagnostics/schemas.py`
- Create: `tests/regression/test_pn2d_minimal6_diagnostic_contracts.py`

**Interfaces:**
- Produces: `QuantityRecord`, `StateIdentity`, `SourceKind`, `SupportKind`, `convert_value`, `classify_pair`, `validate_formula_difference_v1`, and `validate_bv_comparison_v1`.

- [ ] **Step 1: Write RED tests for identities, units, zeros, and schemas**

Cover SI conversions for `V/cm`, `A/cm^2`, `cm^-3`, `cm^2/(V s)`, `cm^-1`, `cm/s`, and native generation units; reject dimensionally invalid conversions. Require unique keys `(run_id, topology, bias_V, carrier, support_kind, support_id, quantity, source, formula_version)`.

- [ ] **Step 2: Define immutable ledger records**

Use enums for support and source identity. Require each record to carry value, unit, sign convention, raw source path/hash, formula version, state identity, and optional geometric-zero reason. Do not use free-form source labels in analysis code.

- [ ] **Step 3: Define `vela.pn2d_minimal6_formula_difference.v1`**

Require input/audit provenance, exact state matrix, row counts, waterfall path definitions, interaction records, dominance rules, named residual, artifact hashes, and the diagnostic-only disclaimer.

- [ ] **Step 4: Define `vela.pn2d_minimal6_bv_comparison.v1`**

Require per-solver/per-topology accepted and failed transitions, exact checkpoint identities, terminal current, maximum field, native/reconstructed source integrals, convergence metadata, curve artifact hashes, and the not-physical-BV disclaimer.

- [ ] **Step 5: Implement deterministic CSV/JSON serialization**

Sort by topology, numeric bias, carrier, support kind, canonical support ID, quantity, source, and formula version. Serialize non-finite values as validation errors, never JSON `NaN`.

- [ ] **Step 6: Run focused tests**

Run `python -m unittest tests.regression.test_pn2d_minimal6_diagnostic_contracts -v`; expected result is PASS.

### Task 3: Independent Physics And Support-Conversion Library

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/geometry.py`
- Create: `scripts/pn2d_minimal6_diagnostics/physics.py`
- Create: `scripts/pn2d_minimal6_diagnostics/support.py`
- Create: `tests/regression/test_pn2d_minimal6_diagnostic_physics.py`

**Interfaces:**
- Produces: `p1_gradient`, `project_vector_to_edge`, `van_overstraeten_alpha`, `invert_alpha`, `infer_ni_eff`, `integrate_nodal_field`, `integrate_cell_field`, `map_local_sources_to_nodes`, and explicit node/edge/cell conversion operators.

- [ ] **Step 1: Write analytic RED tests**

Use affine scalar fields on both canonical topologies to test P1 gradients and orientation invariance. Test signed current projections in both edge directions, Van Overstraeten low/high branches and inverse-domain rejection, OldSlotboom `ni_eff` inference, partial-volume truncation, and conservative node mapping.

- [ ] **Step 2: Implement geometry independently of C++ CSV intermediates**

Construct shape-function gradients from canonical coordinates and triangle tuples. Validate positive signed area and canonicalize IDs only after validating the original CCW tuple.

- [ ] **Step 3: Implement alpha forward/inverse evaluation**

Copy parameter values from the tracked Sentaurus `models.par` and independently check them against Vela's production implementation. Return all valid inverse candidates when the piecewise alpha curve is not one-to-one; choose no driving-field hypothesis silently.

- [ ] **Step 4: Implement `ni_eff` diagnostics**

Infer effective intrinsic density from each internally consistent `(psi, phi_qf, n/p)` relation with documented sign conventions. Emit electron/hole estimates and a consistency residual; do not average inconsistent values into a single authoritative number.

- [ ] **Step 5: Implement named support conversions**

Every node-to-cell, node-vector-to-edge, edge-to-cell, and local-edge-to-node conversion returns both values and normalized weights. Assert conservation for extensive quantities and prohibit averaging native `AvalancheGeneration` before integration.

- [ ] **Step 6: Verify Python and C++ control formulas**

Run the focused Python tests plus `build-release/test_impact_ionization.exe` and `build-release/test_cell_reconstructed_avalanche.exe`. Expected: all pass and the existing `<5e-12` formula gate remains unchanged.

### Task 4: Fixed-State Formula-Difference CLI And Counterfactual Engine

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/counterfactual.py`
- Create: `scripts/diagnose_pn2d_minimal6_formula_difference.py`
- Create: `tests/regression/test_pn2d_minimal6_formula_difference.py`

**Interfaces:**
- CLI: `--state-root`, `--audit-root`, and `--out-dir`.
- Produces: `quantity_ledger.csv`, `factor_waterfall.csv`, `root_cause_summary.json`, and `root_cause_summary.md`.

- [ ] **Step 1: Write RED end-to-end and adversarial tests**

Use the committed synthetic six-state fixture augmented with deterministic native-source/velocity/integral fields. Cover missing fields, wrong units, inexact bias, reversed topology, duplicate nodes, hash mutation, and reconstructed source falsely labelled native.

- [ ] **Step 2: Build the ledger from raw and replayed inputs**

Emit node state/BGN quantities; cell `-grad(psi)`, `grad(phi_n)`, and `grad(phi_p)` components/magnitudes/directions; edge signed current projections, magnitudes, SG flux, midpoint density, mobility, and impact field; exported/recomputed alpha; and all three source families.

- [ ] **Step 3: Define the two source anchors**

Use native integrated `AvalancheGeneration` as the external anchor. Store Sentaurus `alpha*|J|/q` and Vela `alpha*flux*partial_volume` as reconstructed anchors with different `SourceKind` values. The native-minus-reconstructed Sentaurus gap is part of `sentaurus_internal_semantics_residual`.

- [ ] **Step 4: Implement forward and reverse counterfactual paths**

Evaluate factors in the declared order and its reverse: `ni_eff/BGN`, gradient recovery, mobility, current semantics, impact driving field, alpha law, partial volume, and source-to-node mapping. Each step replaces exactly one named operator and re-evaluates downstream dependent quantities.

- [ ] **Step 5: Add 2x2 interaction matrices**

For every adjacent pair whose absolute forward/reverse contribution differs by more than `0.3 dex`, evaluate baseline, A-only, B-only, and A+B. Store the interaction term `AB - A - B + baseline` with sign and path identity.

- [ ] **Step 6: Implement exact closure and dominance scoring**

For every nonzero `-12/-19 V` state, assert `native log gap = named path contributions + residual` within `1e-10 dex`. Score dominance only across both biases and both topologies; never include `0 V` or geometric-zero states.

- [ ] **Step 7: Run focused CLI acceptance**

Require exact `36/54/24` base identities, deterministic output hashes on two runs, schema validation, and explicit residual for every eligible gap.

### Task 5: Formula-Difference Figures And Root-Cause Documentation

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/plots.py`
- Modify: `scripts/diagnose_pn2d_minimal6_formula_difference.py`
- Modify: `tests/regression/test_pn2d_minimal6_formula_difference.py`
- Modify: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`

- [ ] **Step 1: Add figure contract tests**

Require a fixed manifest of gradient, current/alpha, source waterfall, interaction, and topology-symmetry PNG/PDF pairs. Validate signatures, decoding, dimensions, nonblank pixels, units, and disclaimer text.

- [ ] **Step 2: Plot without hiding zeros or direction**

Use signed projections where relevant; use explicit geometric-zero markers instead of log-floor points. Separate native Sentaurus source from both reconstructions visually.

- [ ] **Step 3: Link root causes to implementation entry points**

The Markdown summary must link each classified factor to the tracked Sentaurus deck/parameter entry, independent Python function, and exact C++ production helper or assembler entry. Parameter agreement remains a control unless data disproves it.

- [ ] **Step 4: Perform manual QA**

Inspect node/edge/triangle identities, mirror symmetry, carrier signs, units, zero classification, waterfall closure, and native/reconstructed labelling. Record the inspection checklist and reviewer/date in the manifest.

### Task 6: Minimal6 Vela And Sentaurus Diagnostic-Sweep Drivers

**Files:**
- Create: `scripts/run_pn2d_minimal6_diagnostic_sweep.py`
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_sweep_sdevice.cmd`
- Create: `tests/regression/test_pn2d_minimal6_diagnostic_sweep.py`

**Interfaces:**
- Produces independent Vela and Sentaurus sweep roots plus checkpoint manifests; it does not alter the six-state authoritative schema.

- [ ] **Step 1: Write RED orchestration tests**

Test exact integer target generation `0,-1,...,-20`, separate topology roots, topology-matched mesh/doping, immutable physics configuration, state prefix uniqueness, accepted endpoint validation, failed-transition retention, and refusal to interpolate a missing checkpoint.

- [ ] **Step 2: Generate Vela decks from one immutable template**

Use canonical `mesh.json`, existing `node_doping_file`, fixed-state audit model choices, and `write_state_every_point_prefix`. Deep-compare every generated deck against the template and allow changes only to topology input paths, segment start/end bias, restart state, and output paths.

- [ ] **Step 3: Run segmented Vela continuation**

For each topology, solve `[0,-1]`, then restart sequentially through `[-19,-20]`. Preserve the first failing transition, solver exit code, stdout/stderr, Newton/continuation diagnostics, and the deepest accepted exact checkpoint. Do not relax gates or retry with changed physics.

- [ ] **Step 4: Add the independent Sentaurus sweep deck**

Use the same explicit topology, doping, and tracked model file as fixed-state export. Write a new sweep manifest and exact integer checkpoint TDRs; do not modify or merge the authoritative six-state manifest.

- [ ] **Step 5: Record per-point observables**

For each accepted point record both terminal currents with sign convention, maximum field, native and reconstructed avalanche-source integrals, state hash, and convergence metadata. Record the first rejected transfer as a row with no fabricated physical observables.

- [ ] **Step 6: Run mocked acceptance before expensive solvers**

Use fake runners to exercise full success, Vela early failure, Sentaurus early failure, wrong bias, missing state, and hash tampering. Expected: schemas pass only for internally complete evidence packages.

### Task 7: Sweep Comparison And Fixed-State Root-Cause Recheck

**Files:**
- Create: `scripts/compare_pn2d_minimal6_diagnostic_sweeps.py`
- Create: `tests/regression/test_pn2d_minimal6_sweep_comparison.py`

**Interfaces:**
- Produces: `sweep_comparison.csv`, `sweep_comparison.json`, `sweep_comparison.md`, and I-V/source/field/topology figures conforming to `vela.pn2d_minimal6_bv_comparison.v1`.

- [ ] **Step 1: Write RED comparison tests**

Cover terminal-current sign alignment, one-volt growth ratios, max-field/source comparison, deepest-common-bias calculation, topology sensitivity, missing-tail handling, and failure-transition reporting.

- [ ] **Step 2: Compare only exact common checkpoints**

Never interpolate solver curves. Calculate ratios only when both points are accepted and nonzero; otherwise emit a typed unavailable/zero classification.

- [ ] **Step 3: Re-run the quantity ledger at `0/-12/-19 V`**

Feed each solver's self-consistent exact checkpoint through the same ledger normalization. Re-evaluate the fixed-state factor ranking with the same path and dominance rules and report whether each fixed-state dominant factor remains dominant, changes rank, or becomes unidentifiable.

- [ ] **Step 4: Emit diagnostic-only figures and summary**

Plot terminal currents, 1 V growth ratios, maximum fields, source integrals, and sketch/mirror ratios. Mark solver termination explicitly and never extrapolate a breakdown voltage.

- [ ] **Step 5: Validate schema and closure**

Require every eligible fixed-state and sweep gap to have named contributions plus residual, all artifact hashes to verify, and both solver configurations to be embedded or hash-addressed.

### Task 8: Full Verification And Final Evidence Package

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md`
- Modify: `docs/validation/pn2d_bv_current_progress_summary.md`
- Modify: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`

- [ ] **Step 1: Run all Minimal6 Python regressions**

Run the existing 59-test group plus all new diagnostic contract, physics, formula-difference, sweep, and comparison modules. Record exact counts and elapsed time.

- [ ] **Step 2: Build and run C++ regressions**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target `
  pn2d_minimal6_operator_audit `
  test_fixed_state_operator_audit `
  test_impact_ionization `
  test_cell_reconstructed_avalanche --parallel 2
build-release\test_fixed_state_operator_audit.exe
build-release\test_impact_ionization.exe
build-release\test_cell_reconstructed_avalanche.exe
```

- [ ] **Step 3: Run the real fixed-state report**

Require recovered/regenerated authoritative inputs, successful replay, exact row counts, parity/formula gates, figure QA, schema validation, waterfall closure, interaction records where triggered, and a root-cause ranking or an explicit insufficient-data conclusion.

- [ ] **Step 4: Run both real diagnostic sweeps**

Accept full `-20 V` completion or a preserved first failure from either solver as a valid experimental outcome. Do not weaken configuration to turn a failure into a pass.

- [ ] **Step 5: Audit manifests and documentation**

Verify executable, source, deck, configuration, raw input, state, CSV, JSON, Markdown, and figure hashes. Confirm every claim in validation docs points to a manifest identity and no generated artifact is staged.

- [ ] **Step 6: Run final repository checks**

Run `git diff --check` and `git status --short`. Expected: no whitespace errors; only intended tracked scripts, tests, schemas, decks, and validation documents are changed.

## Stop Conditions

- If neither an authentic archive nor the `sentaurus` target can produce the six exact extended-field states, stop after Task 1 and publish only the data-recovery failure evidence.
- If the installed Sentaurus version cannot export a required native field, stop and revise the field contract explicitly; do not reconstruct it under a native label.
- If replay provenance, state parity, formula parity, topology, units, schema, or hash checks fail, do not generate root-cause rankings.
- If either nonlinear solver terminates early, retain and compare the common accepted prefix; do not tune physics or claim a physical breakdown voltage.

## Commit Boundaries

Commit after each task with narrowly scoped messages: extended state contract, diagnostic core, physics oracle, fixed-state counterfactual report, figure/report QA, sweep orchestration, sweep comparison, and final validation evidence. Generated `build-release/` artifacts remain uncommitted.
