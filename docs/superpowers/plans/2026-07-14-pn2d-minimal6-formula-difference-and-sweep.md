# PN2D Minimal6 Formula-Difference And Diagnostic-Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce reproducible fixed-state causal evidence for the Sentaurus/Vela avalanche-source gap, then test whether the same factors remain dominant in independent sketch/mirror minimal6 nonlinear diagnostic sweeps without changing production physics or solver defaults.

**Architecture:** First recover or regenerate authoritative Sentaurus inputs and fail closed on provenance. Put units, topology, formula reconstruction, counterfactual evaluation, schemas, and plotting behind a small reusable Python package; keep the fixed-state and sweep CLIs as orchestration layers. Treat fixed-state diagnosis and nonlinear sweep as separate publishable phases joined only through versioned ledger records and explicit state identities.

**Tech Stack:** Python 3 standard library plus existing NumPy/Matplotlib/Pillow dependencies, C++20/CMake/Ninja, Catch2, MSYS2 UCRT64, Sentaurus Device 2018 through the existing `sentaurus` SSH workflow.

## Background And Problem Statement

### Background

- The existing Minimal6 fixed-state audit was built to compare Sentaurus and Vela on the same six-node states without involving nonlinear branch selection. Historical handoff evidence reports exact `36 node / 54 edge / 24 triangle` identities, imported-state parity below `1e-12`, and same-formula C++/Python error below `5e-12`; this evidence must be replayed from sealed inputs before it is treated as current.
- That audit establishes topology, state, transport, and formula replay consistency. It does not explain the remaining native avalanche-source magnitude gap, and it does not establish a self-consistent Vela BV branch.
- A local recovery candidate now contains sketch/mirror states at `0/-12/-19 V`, including raw Sentaurus `ImpactIonization`, scalar carrier-speed fields, and ionization-integral fields. Because the directory is ignored and its original manifest is v1, existence alone is not sufficient provenance; Task 1 must seal it as an immutable v2 evidence package or regenerate it under a new run ID.
- Separate coarse7x3 evidence showed that reaching `-20 V` and converging Newton can still leave terminal current at leakage level. Therefore nonlinear convergence, density shape, and endpoint bias are recorded separately from multiplication-current branch classification.

### Questions This Plan Must Answer

1. For exact fixed Sentaurus states, how much of the integrated native-source log gap is attributable to `ni_eff/BGN`, gradient recovery, mobility, current semantics, impact driving field, alpha law, geometric volume, and source-to-node mapping?
2. Do forward and reverse replacement paths agree closely enough to support a dominant-factor statement, and what interactions or residual remain?
3. At exact common nonlinear checkpoints, do the same factors remain dominant, change rank, or become unidentifiable?
4. Did each solver merely converge at the requested bias, or did its terminal current reach the declared multiplication-like order relative to the same-bias Sentaurus reference?

### In Scope

- Reproducible state recovery, schema validation, hashing, independent formula reconstruction, support conversion, causal counterfactuals, diagnostic plotting, Vela/Sentaurus Minimal6 sweeps, and validation documentation.
- Tracked Python tooling, tests, schemas, Sentaurus decks, Vela templates, and Markdown evidence. Generated solver states and reports remain ignored artifacts.

### Out Of Scope

- Changing production C++ physics, solver defaults, calibration parameters, source mapping, or public solver APIs.
- Claiming a physical breakdown voltage, fitting production parameters to Minimal6, interpolating missing checkpoints, or using a reconstructed source under a native-source label.
- Treating fast Newton convergence, a multiplication-like density profile, or arrival at `-20 V` as sufficient branch-recovery evidence.

### Evidence Priority

1. Sealed raw inputs and manifests with verified hashes.
2. Exact-bias replay and schema/parity/formula gates.
3. Deterministic derived ledgers, counterfactual closure, and figures.
4. Historical handoff values only as context; they never override a failed current-run gate.

## Technical Review And Revisions

1. A local recovery candidate now exists at `build-release/reference_tcad/pn2d_sentaurus2018_minimal6/state_exports/minimal6_states_live_20260713_v2`. Its current manifest reports `outputs_complete=true`, six passed states, and SHA-256 `b44ad95d5df6d57383ba3d5b292818568e358d67f0fc0424ee72f95b673e8aaa`. Treat this as a candidate to validate and seal, not as authoritative solely because it exists.
2. Sentaurus deck keyword `AvalancheGeneration` is exported by the installed 2018 workflow under raw field name `ImpactIonization` with unit `cm^-3*s^-1`. The normalized ledger name is `sentaurus_native_avalanche_generation`; Vela's `AvalancheGeneration` remains a distinct raw source. Sentaurus `alpha * |J| / q` is a separately labelled reconstruction.
3. A single ordered waterfall is path-dependent. The report must emit forward and reverse paths, and emit a 2x2 matrix for adjacent factors whose path contributions differ by more than `0.3 dex`.
4. BGN, mobility, current semantics, driving field, alpha law, geometric volume, and node mapping do not all live on the same support. Every substitution must pass through an explicitly named support-conversion operator; direct node/edge/triangle mixing is forbidden.
5. Sentaurus ionization integrals are path-dependent diagnostics, not local alpha or local generation. Store them in the ledger, but do not use them as interchangeable waterfall factors.
6. Do not grow `scripts/audit_pn2d_minimal6_fixed_state.py` into a second analysis framework. Preserve it as the imported-state parity/replay gate and build the new causal analysis from focused modules.
7. Vela's current `node_doping_file` requires `node_id`, `donors_cm3`, and `acceptors_cm3`. The ledger derives net doping from those raw columns; every manifest records raw donor/acceptor hashes, the derived-net formula, and the representation consumed by each solver.
8. Integer checkpoint orchestration must not imply that Vela internally used 1 V steps. The driver may use adaptive continuation inside each `[k, k-1] V` segment, but only a solver-accepted endpoint exactly equal to the integer target may be published.
9. Do not add new required fields to `vela.pn2d_minimal6_states.v1`. Extended fixed states use `vela.pn2d_minimal6_states.v2`; report schemas remain `vela.pn2d_minimal6_formula_difference.v1` and `vela.pn2d_minimal6_bv_comparison.v1`.
10. The recovered Sentaurus 2018 manifests expose `eVelocity` and `hVelocity` as one-component scalar speed fields. The plan accepts those scalars for the first causal pass. A vector-velocity requirement is a separate capability change and must not be inferred from the deck keyword.
11. The work has two execution phases. Phase A (Tasks 1-5) produces a publishable fixed-state report. Phase B (Tasks 6-8) cannot start until the Phase A gate passes and produces a separate diagnostic-sweep package.

## Global Constraints

- Do not modify production C++ APIs, production formulas, default model parameters, source mapping, or solver defaults.
- Do not commit generated Sentaurus, Vela, plot, or report artifacts; keep them below ignored `build-release/` roots.
- Use `D:\msys64\ucrt64\bin` and `D:\msys64\usr\bin` first on `PATH`.
- At execution time, create a fresh branch worktree below ignored `.worktrees/`; do not reuse or remove the dirty `pn2d-source-factor-bisection` worktree.
- Preserve `vela.pn2d_minimal6_states.v1` compatibility. Only `v2` may require the extended native-source, scalar-speed, and ionization-integral fields.
- Required fixed-state matrix is exactly `sketch/mirror x 0/-12/-19 V`, with `6 nodes / 9 unique edges / 4 CCW triangles` per state.
- Required fixed-state row counts remain exactly `36 node / 54 edge / 24 triangle`.
- Imported-state parity must be `<1e-12`; same-formula C++/Python hybrid error must be `<5e-12`.
- Bias identity tolerance is `1e-12 V`; nearest-state and interpolated-state substitution are forbidden.
- Geometry support `<=1e-27 m^2` or both compared source magnitudes `<=1e-285` is classified as `geometric_zero`; it must not receive an artificial log ratio.
- `0 V` is a zero-field/numerical-floor control and is excluded from dominant-root-cause scoring.
- A factor is dominant only if its signed direction agrees at `-12 V` and `-19 V` for both topologies and its symmetric contribution, defined as the mean of forward and reverse path contributions, explains at least 50% of the absolute native-source log gap.
- Every external log gap must equal named factor contributions plus a named `sentaurus_internal_semantics_residual`, within the report serialization tolerance.
- If the absolute residual exceeds 25% of the absolute native-source log gap in any scored state, the report must emit `insufficient_data` instead of a dominant-root-cause claim.
- Every curve and report must say `minimal6 diagnostic sweep; not a physical BV curve`.

## Execution Phases And Gates

**Phase A - fixed-state causal diagnosis (Tasks 1-5):** local-only work may proceed through mocked and synthetic tests. Real reporting requires a sealed `vela.pn2d_minimal6_states.v2` root, exact `36/54/24` identities, parity and formula gates, deterministic report hashes, and exact closure. A residual at or below 25% permits dominance scoring; a larger residual produces valid `insufficient_data` when all integrity gates still pass.

**Phase B - nonlinear diagnostic sweeps (Tasks 6-8):** starts only after Phase A passes or explicitly concludes `insufficient_data` with valid evidence. A solver reaching an integer bias is not proof that it recovered the multiplication-current branch; every accepted point must also carry a typed branch classification based on terminal-current order.

Remote Sentaurus execution is an external gate, not an implicit step. Always try the named local recovery candidate first. If regeneration is required, preserve the failed local validation report and use a new run ID; never overwrite `minimal6_states_live_20260713_v2`.

## Execution Preflight

1. Detect isolation before creating anything: run `git rev-parse --path-format=absolute --git-dir`, `git rev-parse --path-format=absolute --git-common-dir`, `git rev-parse --show-superproject-working-tree`, and `git branch --show-current`. If git-dir differs from common-dir and this is not a submodule, reuse the existing linked worktree.
2. If execution starts from the normal checkout, obtain user consent and prefer the platform's native worktree command. Only when no native worktree command exists, verify `git check-ignore -q .worktrees/probe` and create the fallback with:

```powershell
git worktree add .worktrees\pn2d-minimal6-formula-sweep `
  -b codex-pn2d-minimal6-formula-sweep
```

3. In the isolated worktree, establish a clean baseline before Task 1:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build build-release --target pn2d_minimal6_operator_audit --parallel 2
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_sentaurus_gate `
  tests.regression.test_pn2d_minimal6_fixed_state_audit `
  tests.regression.test_pn2d_minimal6_node_quantity_comparison -v
```

Expected: configure/build exit code 0 and the existing Minimal6 group reports `OK`. If the baseline fails, stop and report the failure before plan implementation.


## Task Contract Summary

The table below is normative. The detailed steps later in the plan explain how to satisfy each contract; when wording differs, the stricter provenance, test, or stop requirement wins.

| Task | Entry condition | Required work | Required deliverables | Task completion condition | Task-local stop condition |
|---|---|---|---|---|---|
| 1. State recovery | Preflight baseline is green; named local candidate is readable or remote regeneration is available. | Preserve v1, validate real raw names/units/components, seal all files and hashes as v2, replay fixed-state audit. | Exporter/test changes, `states.v2` schema, sealed six-state package, validation document. | Six exact sketch/mirror `0/-12/-19 V` states pass topology, bias, field, hash, `36/54/24`, parity, and formula gates. | Stop Phase A if neither local sealing nor a new remote run supplies authentic required fields; never synthesize or relabel them. |
| 2. Ledger and schemas | Task 1 completion condition passes and committed synthetic fixtures are available. | Implement typed identities, units, zero classification, deterministic serialization, and three report/manifest schemas. | Diagnostics core, three JSON schemas, golden and invalid fixtures, focused tests. | Every schema accepts its golden fixture, rejects invalid fixtures, all unit/identity tests pass, and JSON contains no non-finite values. | Stop downstream analysis if identities are non-unique, dimensions conflict, or serialization is nondeterministic. |
| 3. Independent physics | Task 2 types and units are stable; tracked geometry and `models.par` are available. | Implement independent geometry, alpha forward/inverse, `ni_eff`, integration, and named support conversions. | Physics/support modules and analytic tests. | Analytic invariants, conservation checks, Python controls, and both C++ control executables pass; formula error stays below `5e-12`. | Stop counterfactual work on orientation, conservation, parameter-hash, inverse-domain, or formula-control failure. |
| 4. Fixed-state diagnosis | Tasks 1-3 pass and a sealed real or committed synthetic six-state root is supplied. | Build the ledger, source anchors, dependency DAG, forward/reverse paths, interactions, closure, and dominance scoring. | Ledger CSV, waterfall CSV, JSON/Markdown root-cause summaries, focused tests. | Two identical runs have identical hashes; all eligible gaps close within `1e-10 dex`; output is either a valid dominant-factor result or typed `insufficient_data`. | Stop ranking if provenance, parity, units, topology, hashes, or closure fail. A residual above 25% yields `insufficient_data`, not a hard failure. |
| 5. Figures and Phase A evidence | Task 4 outputs validate. | Produce fixed figure manifest, visualize zeros/direction honestly, link factors to source code, perform manual QA. | PNG/PDF pairs, figure manifest, completed validation document. | Figure contract, visual QA, schema, hashes, labels, and Phase A completion gate pass; evidence is committed before Phase B. | Do not enter Phase B for missing provenance, invalid units, parity/schema failure, or unclosed gaps. A typed `insufficient_data` caused by residual size, interactions, or inconsistent factor ranking may proceed. |
| 6. Sweep drivers | Phase A gate permits Phase B. | Implement immutable Vela template, segmented continuation, Sentaurus sweep deck, checkpoint manifests, failure retention, and branch classification. | Vela template, Sentaurus deck, driver, sweep schema usage, mocked tests. | All mocked success/failure/tamper/branch cases pass and generated configurations differ only in allowed fields with verified hashes. | Stop real execution on configuration drift, wrong topology/doping, inexact endpoint identity, schema failure, or fabricated observables. |
| 7. Sweep comparison | Task 6 mocked contracts pass; two valid synthetic or real sweep manifests exist. | Compare exact common checkpoints, re-run ledgers at `0/-12/-19 V`, calculate typed ratios/classes, render diagnostic figures. | Comparison CSV/JSON/Markdown and figure set. | Missing tails and solver failures remain explicit; no interpolation/extrapolation occurs; schemas, hashes, closure, and branch inputs validate. | Stop comparison claims if no exact common checkpoint or if configuration/state identity cannot be verified; still publish the validation failure. |
| 8. Final evidence | Tasks 1-7 implementation tests pass, real inputs and executable solvers are available, and no S1-S4/S6 hard stop is active. | Run full Python/C++ verification, real fixed-state report, real sweeps, manifest audit, and documentation update. | Updated validation/progress documents and hash-addressed evidence package. | All applicable tests/builds pass; Phase A result is valid; Phase B records either full completion or first solver failure plus common-prefix comparison; repository checks are clean. | Use the global terminal-state rules below. A preserved nonlinear early failure is valid experimental evidence, not by itself a plan failure. |
---

### Task 1: Authoritative-State Recovery And Extended Field Contract

**Files:**
- Inspect: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd`
- Modify: `scripts/export_pn2d_minimal6_states.py`
- Modify: `tests/regression/test_pn2d_minimal6_state_export.py`
- Create: `schemas/vela.pn2d_minimal6_states.v2.schema.json`
- Create: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`

**Interfaces:**
- Consumes: the exact local candidate `build-release/reference_tcad/pn2d_sentaurus2018_minimal6/state_exports/minimal6_states_live_20260713_v2`, or a newly generated run from the controlled `sentaurus` SSH target.
- Produces: a sealed `vela.pn2d_minimal6_states.v2` root with six exact states, immutable member hashes, and explicit raw-to-normalized field identities; the v1 contract remains unchanged.

- [ ] **Step 1: Write RED contract tests**

Add named tests for raw `ImpactIonization` as normalized `sentaurus_native_avalanche_generation`, scalar `eVelocity`/`hVelocity`, scalar electron/hole/mean ionization integrals, explicit units, complete global-node mapping, v1 compatibility, exact bias identity, and manifest/member hash mutation. The installed contract is `(ImpactIonization, 1, cm^-3*s^-1)`, `(eVelocity/hVelocity, 1, cm*s^-1)`, and `(eIonIntegral/hIonIntegral/MeanIonIntegral, 1, 1)`.

- [ ] **Step 2: Verify the installed SDevice plot contract**

Confirm that the tracked deck already requests `AvalancheGeneration`, `eVelocity`, `hVelocity`, `eIonIntegral`, `hIonIntegral`, and `MeanIonIntegral`. Do not edit it merely to rename exported fields. If a regenerated run exposes a different raw name or component count, fail the contract and record the capability result; never synthesize a vector or native field.

- [ ] **Step 3: Add the v2 schema and validation without changing v1**

Represent every v2 field as `(normalized_name, raw_name, components, unit, semantic_role)`. Require exactly one region-0 match and reject duplicate aliases. Preserve raw artifact paths and SHA-256 values in each state record, add a root source-manifest hash, and validate `schemas/vela.pn2d_minimal6_states.v2.schema.json` before sealing.

- [ ] **Step 4: Validate the named local candidate before remote execution**

Validate only the named local candidate first. Require source manifest SHA-256 `b44ad95d5df6d57383ba3d5b292818568e358d67f0fc0424ee72f95b673e8aaa`, all member hashes, six-state identity, exact biases, topology, fields, and `outputs_complete=true`. Never rewrite its v1 manifest; write any sealed v2 package under a new sibling run ID.

- [ ] **Step 5: Regenerate under a new run ID if recovery fails**

Set `$runId = "minimal6_states_v2_{0}" -f (Get-Date -Format "yyyyMMdd_HHmmss")`, then run `scripts/export_pn2d_minimal6_states.py --topologies sketch,mirror --biases 0,-12,-19 --run-id $runId --ssh-target sentaurus`. A network, authentication, syntax, or field-capability failure stops Phase A after preserving the partial manifest and logs.

- [ ] **Step 6: Replay the C++ audit and seal provenance**

For every state, generate topology-matched `mesh.json` and immutable `audit.json`, run `build-release/pn2d_minimal6_operator_audit.exe`, and record executable/config/input/output SHA-256 values and the full argv array. Then run the existing fixed-state CLI and require `36/54/24`, parity `<1e-12`, and formula error `<5e-12`.

- [ ] **Step 7: Verify focused regressions**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_fixed_state_audit -v
```

Expected: all tests pass; v1 regressions remain unchanged; the validation document names the sealed v2 run, source v1 hash, exact raw-field identities, and whether local recovery or regeneration supplied the inputs.

### Task 2: Shared Quantity-Ledger Core And Report Schemas

**Files:**
- Create: `scripts/pn2d_minimal6_diagnostics/__init__.py`
- Create: `scripts/pn2d_minimal6_diagnostics/contracts.py`
- Create: `scripts/pn2d_minimal6_diagnostics/units.py`
- Create: `scripts/pn2d_minimal6_diagnostics/ledger.py`
- Create: `scripts/pn2d_minimal6_diagnostics/schemas.py`
- Create: `schemas/vela.pn2d_minimal6_formula_difference.v1.schema.json`
- Create: `schemas/vela.pn2d_minimal6_bv_comparison.v1.schema.json`
- Create: `schemas/vela.pn2d_minimal6_sweep_manifest.v1.schema.json`
- Create: `tests/regression/test_pn2d_minimal6_diagnostic_contracts.py`

**Interfaces:**
- Produces: `QuantityRecord`, `StateIdentity`, `SourceKind`, `SupportKind`, `BranchKind`, `convert_value`, `classify_pair`, `validate_formula_difference_v1`, `validate_bv_comparison_v1`, and `validate_sweep_manifest_v1`.

- [ ] **Step 1: Write RED tests for identities, units, zeros, and schemas**

Cover SI conversions for `V/cm`, `A/cm^2`, `cm^-3`, `cm^2/(V s)`, `cm^-1`, `cm/s`, and native generation units; reject dimensionally invalid conversions. Require unique keys `(run_id, topology, bias_V, carrier, support_kind, support_id, quantity, source, formula_version)`.

- [ ] **Step 2: Define immutable ledger records**

Use enums for support, source, and branch identity. Require each record to carry value, unit, sign convention, raw source path/hash, formula version, state identity, and optional geometric-zero reason. `BranchKind` is exactly `leakage_like`, `multiplication_like`, or `unidentified`. Do not use free-form identity labels in analysis code.

- [ ] **Step 3: Define `vela.pn2d_minimal6_formula_difference.v1`**

Require input/audit provenance, exact state matrix, row counts, waterfall path definitions, interaction records, dominance rules, named residual, artifact hashes, and the diagnostic-only disclaimer.

- [ ] **Step 4: Define comparison and internal sweep-manifest schemas**

Require per-solver/per-topology accepted and failed transitions, exact checkpoint identities, named contact currents in `A/um`, branch classification and threshold provenance, maximum field, native/reconstructed source integrals, convergence metadata, curve artifact hashes, and the not-physical-BV disclaimer. Define the independent internal checkpoint package in `vela.pn2d_minimal6_sweep_manifest.v1` rather than embedding raw states in the comparison report.

- [ ] **Step 5: Implement deterministic CSV/JSON serialization**

Sort by topology, numeric bias, carrier, support kind, canonical support ID, quantity, source, and formula version. Serialize non-finite values as validation errors, never JSON `NaN`.

- [ ] **Step 6: Run focused tests**

Run `D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_diagnostic_contracts -v`; expected result is PASS, all three JSON schemas accept golden fixtures, and each schema rejects one invalid fixture.

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

Parse the required Van Overstraeten values from the tracked Sentaurus `models.par`, record that file's hash, and evaluate them with an independent Python formula. Do not hard-code a second unversioned parameter table. Check parsed values against Vela's production configuration. Return all valid inverse candidates when the piecewise alpha curve is not one-to-one; choose no driving-field hypothesis silently.

- [ ] **Step 4: Implement `ni_eff` diagnostics**

Infer effective intrinsic density from each internally consistent `(psi, phi_qf, n/p)` relation with documented sign conventions. Emit electron/hole estimates and a consistency residual; do not average inconsistent values into a single authoritative number.

- [ ] **Step 5: Implement named support conversions**

Every node-to-cell, node-vector-to-edge, edge-to-cell, and local-edge-to-node conversion returns both values and normalized weights. Assert conservation for extensive quantities and prohibit averaging native Sentaurus `ImpactIonization` or Vela `AvalancheGeneration` before integration.

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

Use integrated raw Sentaurus `ImpactIonization`, normalized as `sentaurus_native_avalanche_generation`, as the external anchor. Store Sentaurus `alpha*|J|/q` and Vela `alpha*flux*partial_volume` as reconstructed anchors with different `SourceKind` values. Never label either reconstruction as native. The native-minus-reconstructed Sentaurus gap is part of `sentaurus_internal_semantics_residual`.

- [ ] **Step 4: Implement forward and reverse counterfactual paths**

Encode an explicit dependency DAG, then evaluate factors in the declared order and its reverse: `ni_eff/BGN`, gradient recovery, mobility, current semantics, impact driving field, alpha law, partial volume, and source-to-node mapping. Each step replaces exactly one named operator and re-evaluates only its declared downstream dependents. Reject undeclared dependencies.

- [ ] **Step 5: Add 2x2 interaction matrices**

For every adjacent pair whose absolute forward/reverse contribution differs by more than `0.3 dex`, evaluate baseline, A-only, B-only, and A+B. Store the interaction term `AB - A - B + baseline` with sign and path identity.

- [ ] **Step 6: Implement exact closure and dominance scoring**

For every nonzero `-12/-19 V` state, assert `native log gap = named path contributions + residual` within `1e-10 dex`. Define each factor's symmetric contribution as the arithmetic mean of its forward and reverse contributions. Score dominance only across both biases and both topologies; never include `0 V` or geometric-zero states. If any scored state's residual exceeds 25% of its absolute gap, emit `insufficient_data` and no dominant factor.

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

**Phase A completion gate:** run the complete Task 1-5 focused test set, validate the sealed v2 state root and formula-difference schema, and record either a root-cause result or `insufficient_data`. Commit that evidence before starting Task 6. A typed `insufficient_data` result may proceed when provenance, units, topology, parity, formula, schema, deterministic-hash, and closure gates all pass; missing or failed integrity evidence may not proceed.

---

### Task 6: Minimal6 Vela And Sentaurus Diagnostic-Sweep Drivers

**Files:**
- Create: `scripts/run_pn2d_minimal6_diagnostic_sweep.py`
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_sweep_sdevice.cmd`
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/vela/pn2d_minimal6_sweep_template.json`
- Create: `tests/regression/test_pn2d_minimal6_diagnostic_sweep.py`

**Interfaces:**
- Produces independent Vela and Sentaurus sweep roots conforming to `vela.pn2d_minimal6_sweep_manifest.v1`; it does not alter or embed the six-state authoritative schema.

- [ ] **Step 1: Write RED orchestration tests**

Test exact integer target generation `0,-1,...,-20`, separate topology roots, topology-matched mesh/doping, immutable physics configuration, state prefix uniqueness, accepted endpoint validation, failed-transition retention, refusal to interpolate a missing checkpoint, and all three `BranchKind` outcomes.

- [ ] **Step 2: Generate Vela decks from one immutable template**

Create one tracked immutable template using canonical `mesh.json`, a `node_doping_file` with `node_id/donors_cm3/acceptors_cm3`, fixed-state audit model choices, and `write_state_every_point_prefix`. Deep-compare every generated deck against the template and allow changes only to topology input paths, segment start/end bias, restart state, and output paths. Store template and substituted-deck hashes.

- [ ] **Step 3: Run segmented Vela continuation**

For each topology, solve `[0,-1]`, then restart sequentially through `[-19,-20]`. Preserve the first failing transition, solver exit code, stdout/stderr, Newton/continuation diagnostics, and the deepest accepted exact checkpoint. Do not relax gates or retry with changed physics.

- [ ] **Step 4: Add the independent Sentaurus sweep deck**

Use the same explicit topology, doping, and tracked model file as fixed-state export. Write a new sweep manifest and exact integer checkpoint TDRs; do not modify or merge the authoritative six-state manifest.

- [ ] **Step 5: Record per-point observables**

For each accepted point record `Anode` and `Cathode` terminal currents in signed `A/um`, the current extraction formula, maximum field, native and reconstructed avalanche-source integrals, state hash, and convergence metadata. At a common non-geometric-zero checkpoint, classify Vela as `multiplication_like` when `0.1 <= abs(I_vela/I_sentaurus) <= 10`, `leakage_like` when the ratio is `<=1e-3`, and `unidentified` otherwise; store the ratio and threshold version. Record the first rejected transfer as a row with no fabricated physical observables.

- [ ] **Step 6: Run mocked acceptance before expensive solvers**

Use fake runners to exercise full success, Vela early failure, Sentaurus early failure, wrong bias, missing state, hash tampering, all branch classifications, and a zero/geometric-zero reference that must produce `unidentified`. Expected: `vela.pn2d_minimal6_sweep_manifest.v1` passes only for internally complete evidence packages.

### Task 7: Sweep Comparison And Fixed-State Root-Cause Recheck

**Files:**
- Create: `scripts/compare_pn2d_minimal6_diagnostic_sweeps.py`
- Create: `tests/regression/test_pn2d_minimal6_sweep_comparison.py`

**Interfaces:**
- Produces: `sweep_comparison.csv`, `sweep_comparison.json`, `sweep_comparison.md`, and I-V/source/field/topology figures conforming to `vela.pn2d_minimal6_bv_comparison.v1`.

- [ ] **Step 1: Write RED comparison tests**

Cover `Anode`/`Cathode` sign alignment, `A/um` normalization, one-volt growth ratios, max-field/source comparison, deepest-common-bias calculation, all branch classifications, topology sensitivity, missing-tail handling, and failure-transition reporting.

- [ ] **Step 2: Compare only exact common checkpoints**

Never interpolate solver curves. Calculate ratios and branch class only when both points are accepted, nonzero, and non-geometric-zero; otherwise emit a typed unavailable/zero classification. Do not relabel an exact-bias `leakage_like` point as branch recovery merely because Newton converged.

- [ ] **Step 3: Re-run the quantity ledger at `0/-12/-19 V`**

Feed each solver's self-consistent exact checkpoint through the same ledger normalization. Re-evaluate the fixed-state factor ranking with the same path and dominance rules and report whether each fixed-state dominant factor remains dominant, changes rank, or becomes unidentifiable.

- [ ] **Step 4: Emit diagnostic-only figures and summary**

Plot terminal currents, 1 V growth ratios, maximum fields, source integrals, and sketch/mirror ratios. Mark solver termination explicitly and never extrapolate a breakdown voltage.

- [ ] **Step 5: Validate schema and closure**

Require every eligible fixed-state and sweep gap to have named contributions plus residual, all artifact hashes to verify, both solver configurations to be embedded or hash-addressed, and every common checkpoint to carry the serialized branch-ratio inputs and threshold version.

### Task 8: Full Verification And Final Evidence Package

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md`
- Modify: `docs/validation/pn2d_bv_current_progress_summary.md`
- Modify: `docs/validation/pn2d_minimal6_formula_difference_2026-07-14.md`

- [ ] **Step 1: Run all Minimal6 Python regressions**

Run the existing Minimal6 group plus all new modules with an explicit command; do not retain a stale expected count after adding tests:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_sentaurus_gate `
  tests.regression.test_pn2d_minimal6_fixed_state_audit `
  tests.regression.test_pn2d_minimal6_node_quantity_comparison `
  tests.regression.test_pn2d_minimal6_diagnostic_contracts `
  tests.regression.test_pn2d_minimal6_diagnostic_physics `
  tests.regression.test_pn2d_minimal6_formula_difference `
  tests.regression.test_pn2d_minimal6_diagnostic_sweep `
  tests.regression.test_pn2d_minimal6_sweep_comparison -v
```

Expected: exit code 0 and `OK`. Record the discovered test count and elapsed time in the validation document.

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

Accept full `-20 V` completion or a preserved first failure from either solver as a valid experimental outcome. Separately report the deepest exact checkpoint classified `multiplication_like`; reaching `-20 V` while `leakage_like` is not branch recovery. Do not weaken configuration to turn a failure into a pass.

- [ ] **Step 5: Audit manifests and documentation**

Verify executable, source, deck, configuration, raw input, state, CSV, JSON, Markdown, and figure hashes. Confirm every claim in validation docs points to a manifest identity and no generated artifact is staged.

- [ ] **Step 6: Run final repository checks**

Run `git diff --check` and `git status --short`. Expected: no whitespace errors; only intended tracked scripts, tests, schemas, decks, and validation documents are changed.


## Completion Conditions

### `COMPLETE`

The implementation is complete only when all of the following are true:

- C1. Tasks 1-8 satisfy every task-completion condition in the contract table, and every tracked change belongs to the declared file set.
- C2. A sealed `vela.pn2d_minimal6_states.v2` root verifies all raw/member hashes, six exact state identities, required field semantics, topology, `36/54/24` row counts, parity `<1e-12`, and same-formula error `<5e-12`.
- C3. Phase A produces schema-valid, deterministic ledgers/reports/figures with closure within `1e-10 dex` and either a supported dominant-factor result or a typed inconclusive result.
- C4. Phase B preserves exact integer checkpoint identity and immutable solver configuration; both real solver attempts produce either a complete `0..-20 V` sweep or an immutable first-failure record. When an exact common prefix exists, comparison artifacts cover that entire prefix without interpolation.
- C5. The full listed Python suite, required C++ build targets/executables, schema checks, hash audit, `git diff --check`, and documentation claim audit pass with current-run evidence.
- C6. Generated states, logs, CSV/JSON reports, and figures remain under ignored `build-release/` roots; only intended source, tests, schemas, decks/templates, and validation documents are tracked.

### `COMPLETE_WITH_INCONCLUSIVE_RESULT`

Scientific inconclusiveness is an allowed completed outcome when implementation and evidence integrity are complete. It applies only to these cases:

- Phase A closes numerically but the semantics residual exceeds 25%, so the result is `insufficient_data` and no dominant factor is named.
- Phase B has valid exact checkpoints but fixed-state factors change rank or become `unidentified` on nonlinear states.
- A nonlinear solver stops early after producing an immutable first-failure record and at least one exact common checkpoint remains available for the comparison required by Task 7.

These cases must not be rewritten as success of a physical model or recovery of a multiplication branch; they are complete diagnostic outcomes with an inconclusive scientific conclusion.

### `NOT_COMPLETE`

Any S1-S4 or S6 condition below makes the plan not complete. S5 is complete only when its first-failure evidence is immutable and Task 7 still has at least one exact common checkpoint; otherwise it is also not complete. Every not-complete run ends as `STOPPED_WITH_EVIDENCE`, not as a scientific conclusion.

## Stop Conditions

### S1. Authoritative-input hard stop

- Trigger: the named local v1 candidate fails hash/identity/field validation and a new remote run cannot produce six exact v2 states.
- Required response: stop after Task 1; preserve the candidate validation report, partial new manifest, exact failing command, stdout/stderr, and remote/local run IDs.
- Forbidden downstream claims: no current fixed-state report, root-cause ranking, Phase B sweep, or branch conclusion.

### S2. Sentaurus field-capability hard stop

- Trigger: Sentaurus 2018 cannot export raw `ImpactIonization`, scalar speeds, or required ionization integrals with the declared unit/component/mapping contract.
- Required response: record the raw names actually observed and revise the schema/plan through review. Do not synthesize a missing field, infer a vector from a scalar, or regenerate under a native label.

### S3. Evidence-integrity hard stop

- Trigger: provenance, state parity, formula parity, topology, units, support conservation, schema, deterministic hash, or waterfall closure fails its stated tolerance.
- Required response: stop the affected task and every dependent task; publish the validation failure and preserve inputs/outputs. Root-cause ranking and comparison claims are prohibited until the failing gate is repaired and rerun.

### S4. Phase-transition hard stop

- Trigger: Phase A ends because of missing provenance, invalid units, topology/parity/formula/schema/hash failure, or unclosed gaps rather than a schema-valid `insufficient_data` result with every integrity gate passed.
- Required response: do not start Task 6. Phase B is allowed only after a fresh Phase A completion gate passes.

### S5. Nonlinear experimental stop

- Trigger: Vela or Sentaurus rejects a transition, exits nonzero, produces an inexact endpoint, or reaches `-20 V` while classified `leakage_like`.
- Required response: preserve the first failing transition and deepest valid exact checkpoint. Never retry with changed physics, relaxed convergence gates, or altered topology/doping under the same run identity.
- Completion effect: early termination is valid experimental evidence only when the manifest is complete and Task 7 can compare at least one exact common checkpoint. Endpoint convergence with `leakage_like` is explicitly branch non-recovery.

### S6. No-common-checkpoint comparison stop

- Trigger: the two solvers have no accepted, identity-verified exact checkpoint in common, or a required state/configuration hash cannot be verified.
- Required response: do not calculate ratios, factor persistence, or branch comparison. Publish manifest-level failure evidence; the final state is `STOPPED_WITH_EVIDENCE`, not complete.

### Universal prohibitions after any stop

- Do not interpolate or substitute a nearby state.
- Do not weaken physics/configuration to convert a stopped run into a pass.
- Do not claim a physical breakdown voltage or multiplication-branch recovery without the explicit current-order gate.

## Commit Boundaries

Commit after each task with narrowly scoped messages: extended state contract, diagnostic core, physics oracle, fixed-state counterfactual report, figure/report QA, sweep orchestration, sweep comparison, and final validation evidence. Generated `build-release/` artifacts remain uncommitted.
