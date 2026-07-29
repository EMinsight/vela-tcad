# PN2D BV process observability and causal-parity development plan

Date: 2026-07-29

Status: WP0-WP2.1 and scientific Task 2 executed; WP3 is next.

Starting point: branch `codex-pn2d-minimal6-operator-audit`, commit
`1350d11`.

Companion scientific acceptance plan:
`docs/superpowers/plans/2026-07-29-pn2d-avalanche-on-bv-knee-parity.md`.

## 1. Purpose and relationship to the parity plan

The companion parity plan defines the scientific target, bias lattices,
quantitative curve gates, stop conditions, mesh requirements, and default-change
policy. This document defines the code and data work required to make those
decisions reproducible.

The immediate engineering objective is not to change the avalanche model. It
is to construct one fail-closed observation pipeline that can identify the
first differing stage in:

```text
psi/QFP -> n,p -> E/grad(QFP) -> mobility/current -> alpha
        -> carrier-split Gava -> geometric integration -> residual/Jacobian
        -> Newton update -> terminal current
```

The pipeline must distinguish:

1. a fixed-state operator difference;
2. a source-to-state nonlinear feedback difference;
3. a continuation or solver-path difference; and
4. an output/support/normalization mismatch.

No work package in this plan authorizes an empirical current, voltage, field,
mobility, lifetime, source, or avalanche-coefficient fit.

## 2. Frozen starting evidence

The following facts are inputs, not hypotheses to be re-tuned:

- M0 has 27 physical nodes and 32 triangles with mesh hash
  `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e`.
- The M0 avalanche-off 21-point Vela/Sentaurus log-current RMSE is
  `6.997e-5 dex`; maximum error is `1.255e-4 dex`.
- Sentaurus implicit avalanche and explicit `GradQuasiFermi` are identical on
  the accepted high-bias lattice.
- The accepted exact high-bias lattice is
  `-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`.
- At the fixed `-20 V` hotspot, the first material change from `-19 V` to
  `-20 V` is electron density, not electric field, QFP-gradient magnitude, or
  mobility.
- The current nonlinear blocker is `nonsmooth_branch_derivative`; the best
  accepted nonzero source-Jacobian analytic/central-FD difference is
  `1.839007985e-6`, above the frozen `1e-8` gate.
- The latest source-scaling and avalanche-off changes do not yet have a sealed
  current-code avalanche-on curve and knee scorecard.

Production choices remain frozen unless a later work package explicitly
authorizes an opt-in candidate:

- BV mobility doping basis: `net_doping`;
- forward-IV mobility doping basis:
  `cell_reconstructed_total_impurity`;
- avalanche model: `van_overstraeten`;
- avalanche drive: `quasi_fermi_gradient`;
- production current approximation: `cell_reconstructed`;
- production source mapping: `triangle_gss_gradqf_truncated`;
- SRH, Old Slotboom BGN, and QFP-gradient high-field mobility: enabled;
- M0 mesh through causal localization; and
- Sentaurus release: `O-2018.06-SP2`.

### 2.1 Execution snapshot after WP0-WP2

The first execution slice completed on 2026-07-29 without changing the Vela
physics or solver:

- WP0 outcome: `implementation_baseline_sealed`. The avalanche-off 21-point
  log-current RMSE is `6.996509075786512e-5 dex`, maximum error is
  `1.254688929143555e-4 dex`, electron/hole source-contact closure is
  `7.48198732275877e-7`/`3.561187670708268e-9`, and terminal-pair closure is
  `3.7621205605551904e-23 A/um`.
- The duplicate avalanche-on Vela runs are deterministic but incomplete. Both
  contain 2,232 accepted states, end at `-19.692187499999644 V`, and fail the
  request for `-19.693749999999643 V` with `max_iterations`. Their normalized
  output hash is
  `dacf644b49080091b09c066bec66b6ccf899afc7671b257fa41ab6968c83d87b`.
  Therefore the baseline is sealed with a deterministic solver-first-failure,
  not with complete `-20 V` curve coverage.
- WP1 outcome: `process_contract_verified`. Thirteen focused tests verify exact
  bias matching, support/provenance/unit declarations, coordinates and
  connectivity, hashes, and fail-closed rejection of unsupported native-edge
  claims and implicit zero fill.
- WP2 outcome: `sentaurus_process_matrix_available`. Two independent VM roots
  each produced 32 exact TDR snapshots, 48,736 field records, and 288
  aggregates on the eight-point high-bias lattice. Input and normalized output
  hashes match, and requested/actual bias error is zero.
- Sentaurus off versus IIC has zero state difference over 7,768 records while
  IIC generation is nonzero (maximum `1.1171098092257102e15 cm^-3 s^-1`).
  Avalanche-on versus `AvalDerivatives` also has zero state difference over
  7,768 records. Native CurrentPlot and Tcl/`ReadMeasure` source integration
  close to `3.977747984452422e-15` relative.
- WP2.1 outcome: `exact_reference_lattices_complete`. Two independent VM
  roots each produced 116 exact snapshots on the 29-point union of the global
  and knee lattices, 176,668 field records, and 1,044 aggregates. Both
  acceptance lattices have no missing rows, requested/actual bias error is
  zero, and normalized/input hashes match.
- WP2.1 off versus IIC and on versus `AvalDerivatives` state differences are
  zero over 28,159 records per comparison. All 348 native/replayed source
  integrals close, with maximum relative error
  `4.830623899888963e-15`.
- Scientific Task 2 outcome: `curve_knee_contract_verified` in 12 synthetic
  scenarios. Current-data classification is `solver_first_failure`; Vela is
  missing global `-20 V` and knee
  `-19.7, -19.8, -19.85, -19.9, -19.95, -20 V`, with the sealed first failure
  at `-19.693749999999643 V`.
- Thirty-six focused Python tests and the complete Release suite
  (`495/495`) pass.

Primary evidence bundles:

- `build-release/pn2d-wp0-implementation-baseline-20260729/acceptance.json`;
- `build-release/pn2d-wp2-process-matrix-pair-20260729/acceptance.json`.
- `build-release/pn2d-wp21-full-lattice-pair-20260729/acceptance.json`.
- `build-release/pn2d-task2-curve-knee-20260730/acceptance.json`.

The source/schema/plan files implementing WP1-WP2 are currently untracked.
They must be reviewed and committed in the boundaries defined in section 15
before WP3 changes are mixed into the worktree. `tmp/` remains user-owned and
must not be staged.

## 3. Architecture and artifact contract

### 3.1 Three paired physics branches

Every target bias must be observable on these branches:

| Branch | Carrier equations | Avalanche evaluation | Avalanche feedback |
|---|---|---|---|
| `avalanche_off` | Poisson + Electron + Hole | disabled | disabled |
| `iic_postprocess` | Poisson + Electron + Hole | enabled with GradQF | excluded from continuity RHS |
| `avalanche_on` | Poisson + Electron + Hole | enabled with GradQF | self-consistent |

The Sentaurus `iic_postprocess` branch is the tutorial's ionization-integrals
with carriers method using `ComputeIonizationIntegrals` and
`AvalPostProcessing`. It must retain `GradQuasiFermi`, not copy the tutorial's
example `EParallel` selector.

The two causal differences are then:

```text
avalanche_off -> iic_postprocess : fixed-state avalanche evaluation
iic_postprocess -> avalanche_on  : nonlinear avalanche feedback
```

ABA/Poisson-only may be added later as a labelled field-path diagnostic. It is
not an acceptance branch because it uses `ElectricField` and is only a
qualitative off-state approximation.

### 3.2 Exact bias sets

Use the companion plan's two lattices:

- global: `0, -1, -2, ..., -20 V`;
- knee: `-18, -18.5, -19, -19.25, -19.5, -19.7, -19.8, -19.85,
  -19.9, -19.95, -20 V`.

The completed WP2 process-pipeline proof used the eight-point accepted
high-bias lattice. It is sufficient to validate the extraction mechanism, but
not to run the final curve/knee gates. WP2.1 must add the missing knee points
`-18.5, -19.25, -19.7, -19.85 V` and the complete global lattice before a
scientific parity result can be emitted. A point must carry its actual contact
voltage from the Current/TDR output. Filename index, requested voltage, minimum
QFP, interpolation, or nearest-state matching is not sufficient evidence.

### 3.3 Per-run directory

Each run writes below an ignored `build-release/` root:

```text
<run-root>/
  manifest.json
  branch_manifest.json
  bundle/
  raw/
  normalized/
    bias_manifest.json
    aggregate.csv
    nodes.csv
    cells.csv
    element_edges.csv
    element_vertices.csv
    contacts.csv
    newton_attempts.csv
    newton_iterations.csv
  analysis/
```

Every normalized row includes:

- simulator, release, branch, requested bias, and actual bias;
- mesh/doping/material/contact/config hashes;
- support kind and stable support identity;
- native, reconstructed, solver-used, or postprocessed provenance;
- carrier, quantity, components, and canonical unit;
- coordinate or connectivity required to remap the support; and
- source file, dataset, and original row/index.

Canonical comparison units are:

- geometry: `um`;
- potential/QFP: `V`;
- density and charge density: `cm^-3`;
- electric field and QFP gradient: `V/cm`;
- mobility: `cm^2/(V s)`;
- velocity: `cm/s`;
- current density: `A/cm^2`;
- avalanche alpha: `cm^-1`;
- generation/recombination: `cm^-3 s^-1`; and
- integrated 2-D current/source: `A/um`.

### 3.4 Delivery roadmap

| Milestone | Work packages | Deliverable | Entry gate to the next milestone |
|---|---|---|---|
| M0: frozen evidence | WP0 | current-code on/off baseline and reproducibility bundle | baseline and duplicate-run gates pass |
| M1: common observations | WP1-WP2.1 | fail-closed schema plus Sentaurus off/IIC/on exact process matrix on both acceptance lattices | every exact state and native aggregate validates; the curve/knee verifier fails closed |
| M2: Vela parity instrumentation | WP3-WP5 | postprocess-only branch, solver-used records, complete nonlinear trace | diagnostic formulas close and do not alter states |
| M3: nonlinear authorization | WP6 | source-only derivative classification and correction, if proven | frozen Jacobian gate passes |
| M4: causal localization | WP7 | first-departure report on two adjacent knee biases | one causal stage has cross-bias evidence |
| M5: correction qualification | WP8 | minimal opt-in candidate, M0/finer-mesh/full-regression evidence | all companion-plan gates and reviews pass |

Priority classification:

- P0: WP0-WP7. These are required to make a defensible causal claim.
- P1: the minimal correction and M0 qualification portion of WP8.
- P2: dose-preserving mesh refinement and any separate production-default
  proposal after the opt-in candidate passes.

## 4. Work package 0 - seal the implementation baseline

### Goal

Prove that later output changes are observational and do not silently alter
the accepted M0 physics.

### Actions

1. Verify `fa1c343` and `1350d11` are ancestors of `HEAD`.
2. Record `git status --short --branch`; preserve the existing untracked
   parity plan and `tmp/` without staging them.
3. Re-run or verify the M0 avalanche-off 21-point baseline and source/contact
   closure.
4. Hash the Vela executable, rendered configuration, mesh, doping, materials,
   and accepted Sentaurus inputs.
5. Capture the present avalanche-on global and knee sweeps twice before any
   code change. If either cannot complete, preserve the same first failed
   transition from both attempts.

### Acceptance

- Avalanche-off RMSE `<= 0.001 dex`, maximum `<= 0.002 dex`.
- Electron/hole source-contact closure `<= 1e-5` when above the source floor.
- Total terminal-pair closure `<= 1e-20 A/um` absolute.
- Duplicate avalanche-on runs have identical accepted bias sequences or the
  same typed first failure.

### Outcome

- `implementation_baseline_sealed`; or
- `baseline_or_determinism_mismatch` and stop.

Executed outcome: `implementation_baseline_sealed`, qualified by
`deterministic_solver_first_failure` at the transition
`-19.692187499999644 -> -19.693749999999643 V`. This qualification is a
coverage fact and must remain visible in every downstream scorecard.

## 5. Work package 1 - define the common process-data schema

### Files

Create:

- `schemas/vela.pn2d_bv_process_run.v1.schema.json`;
- `scripts/pn2d_bv_process_contract.py`;
- `tests/regression/test_pn2d_bv_process_contract.py`.

Modify only if reusable helpers are needed:

- `scripts/pn2d_sentaurus_process_run_contract.py`.

### Actions

1. Define run, branch, bias, support, provenance, field, aggregate, and Newton
   attempt records.
2. Use explicit support enums:
   `physical_node`, `contact_support_vertex`, `cell`, `element_local_edge`,
   `element_local_vertex`, and `contact`.
3. Use explicit provenance enums:
   `native`, `operator_replay`, `reconstructed`, `solver_used`, and
   `postprocessed`.
4. Require `requested_bias_V` and `actual_bias_V`; fail if their difference
   exceeds `1e-10 V` for an exact-state gate.
5. Require file hashes and normalized-output hashes.
6. Reject implicit zero-fill, unknown units, duplicate support keys,
   nonfinite values, missing carrier labels, and nearest-bias substitution.

### Tests

Cover:

- a minimal valid paired run;
- contact-support duplicate vertices separated from physical nodes;
- duplicate and missing bias rows;
- wrong centering or unit;
- unsupported native-edge-current claims;
- nonfinite values;
- hash drift; and
- reordered rows producing the same normalized hash.

### Acceptance

All negative fixtures fail with a typed reason. No consumer needs to infer
support, provenance, unit, or actual bias from a filename.

### Outcome

- `process_contract_verified`; or
- `process_contract_invalid` and stop.

Executed outcome: `process_contract_verified` with 13 focused tests.

## 6. Work package 2 - add exact Sentaurus IIC and spatial snapshots

### Files

Prefer creating a focused orchestrator while reusing existing safe SSH/SCP
helpers:

- create `scripts/run_pn2d_bv_process_matrix_vm.py`;
- modify `scripts/run_pn2d_high_bias_oracle_variant_vm.py` only to extract
  reusable deck-generation helpers;
- reuse `scripts/run_sentaurus_vm_reference.py` for Windows OpenSSH binaries,
  command execution, remote-root validation, and manifests;
- create `tests/regression/test_pn2d_bv_process_matrix_vm.py`.

### Sentaurus deck requirements

Add controlled branches:

1. `avalanche_off`;
2. `iic_postprocess` with:

   ```text
   Math {
     ComputeIonizationIntegrals
     AvalPostProcessing
   }
   ```

3. `avalanche_on`;
4. `avalanche_on_aval_derivatives`, differing only by
   `Math { AvalDerivatives }`, as a numerical-path control.

For every exact target, use an explicit `Quasistationary Goal` segment. At the
end of each segment write:

- a unique non-loadable Plot TDR;
- an optional Save state for exact restart;
- a CurrentPlot row; and
- the runtime Tcl process record.

The Plot list must include, when supported:

- Potential, e/hQuasiFermi, e/hDensity;
- ElectricField and e/hGradQuasiFermi on vertex and element support;
- e/hMobility, e/hVelocity;
- e/h/total current-density vectors;
- e/hAlphaAvalanche;
- e/h/total AvalancheGeneration;
- e/h/MeanIonIntegral for the IIC branch;
- SRHRecombination, Doping, donor/acceptor concentration, and SpaceCharge.

The top-level CurrentPlot section must include:

- electron, hole, and total avalanche-generation integrals;
- maximum ElectricField with coordinates;
- maximum impact ionization with coordinates; and
- the existing Tcl channel for `ReadCoefficient`, `ReadMeasure`, element
  vectors, and carrier-split source reintegration.

Set `CurrentPlot(IntegrationUnit=um Digits=12)` and retain double round-trip
precision in fetched/exported artifacts.

For 2-D `IntegrationUnit=um` generation integrals, canonical conversion to
`A/um` is `q * 1e-12` times the integrated `cm^-3 s^-1` value. The converter
and its tests must keep this factor explicit; no inferred device-depth or
unitless scaling is allowed.

VM execution must tolerate a queued license without duplicating completed
work. The runner shall keep remote completion markers, support `--resume`, and
fetch completed outputs in batches. A client-side SSH timeout alone is not a
typed simulation failure; the remote completion marker and simulator logs
decide the run state.

### Failure-state diagnostics

For the self-consistent branches:

- enable `CNormPrint` for the knee run;
- configure a unique NewtonPlot prefix;
- produce NewtonPlot TDR only for the first failed/rejected transition or a
  predeclared knee probe, not for every successful step; and
- retain complete `.log`, `.plt`, `.tdr`, `.out`, deck, parameter, and Tcl
  files.

### Tests

The dry-run/deck tests must prove:

- IIC differs from avalanche-on only by the declared feedback controls;
- every target is an explicit Goal and has a unique snapshot name;
- GradQF remains selected;
- exact bias and CurrentPlot integration-unit controls are present;
- no ElectricField avalanche substitution is introduced; and
- remote commands are argv-based and use a validated absolute POSIX root.

### Acceptance

- Two independent roots produce identical normalized CurrentPlot and process
  records.
- Every exact point has one TDR, one CurrentPlot row, and one process record.
- IIC generation is nonzero where expected while its carrier state follows
  the avalanche-off branch within solver tolerance.
- `AvalDerivatives` on/off converged states agree within `1e-10` relative, or
  the difference is typed as a solver/continuation branch difference.

### Outcome

- `sentaurus_process_matrix_available`;
- `iic_state_not_decoupled`;
- `exact_snapshot_mismatch`; or
- `sentaurus_solver_path_difference`.

Only the first outcome enters work package 3.

Executed outcome: `sentaurus_process_matrix_available` on the eight-point
high-bias extraction lattice. This verifies WP2 mechanics, but WP2.1 remains
mandatory before scientific curve/knee acceptance.

## 6.1 Work package 2.1 - complete exact acceptance lattices

### Goal

Promote the validated WP2 extraction mechanism from an eight-point engineering
proof to the exact global and knee reference required by the companion plan.

### Actions

1. Add exact Sentaurus off/IIC/on snapshots for the missing knee points
   `-18.5, -19.25, -19.7, -19.85 V`.
2. Add the complete exact global lattice `0, -1, ..., -20 V`, including an
   explicitly labelled equilibrium `0 V` state rather than inferring it from
   the first sweep record.
3. Re-run only missing states using the WP2 resume/completion-marker contract,
   then normalize them into one manifest per branch and lattice.
4. Reproduce new rows in two independent VM roots and compare normalized
   hashes, actual bias, state fields, aggregates, and source reintegration.
5. Connect the completed manifests to scientific Task 2's curve/knee verifier.
   Until all exact rows exist, that verifier must emit
   `incomplete_exact_lattice`, never partial knee metrics.

### Acceptance

- Every global and knee target has exactly one converged exact state for every
  required Sentaurus branch.
- Requested/actual bias difference is `<= 1e-10 V`.
- Duplicate roots have identical normalized hashes and pass the WP2
  state/source checks.
- No interpolation, nearest-state substitution, or implicit zero fill appears
  in an acceptance input.

### Outcome

- `exact_reference_lattices_complete`; or
- `incomplete_exact_lattice` and stop before scientific parity scoring.

Executed outcome: `exact_reference_lattices_complete` on the 29-point union
lattice in two independent VM roots. Scientific Task 2 is authorized.

## 7. Work package 3 - implement a Vela IIC-equivalent observation mode

### Design

Add an opt-in avalanche coupling selector under the existing impact-ionization
configuration:

```json
"impact_ionization": {
  "model": "van_overstraeten",
  "coupling_mode": "self_consistent"
}
```

Allowed values:

- `self_consistent` - current behavior and default;
- `postprocess_only` - compute avalanche process quantities but exclude the
  source from the continuity residual and Jacobian.

Continue to use `model: "none"` for a truly disabled avalanche branch.

### Files

Modify:

- `include/vela/physics/ImpactIonizationModel.h`;
- `src/solver/NewtonSolver.cpp`;
- `src/solver/GummelSolver.cpp`;
- `src/equation/CoupledDDAssembler.cpp`;
- `include/vela/equation/CoupledDDAssembler.h`;
- `src/simulation/DCSweep.cpp`;
- `configs/schema/vela-simulation.schema.json`;
- `docs/config_schema.md`;
- `configs/templates/pn2d_bv.template.json` only if an explicit template
  parameter is needed.

Tests:

- `tests/test_newton_solver.cpp`;
- `tests/test_impact_ionization.cpp`;
- `tests/test_dc_sweep.cpp`.

### Implementation constraints

1. Default-absent and explicit `self_consistent` results must be identical.
2. `postprocess_only` residual and Jacobian must be identical to `model:none`
   when all other physics is the same.
3. `postprocess_only` output must still contain alpha, drive, current, and
   carrier-split generation computed with the declared avalanche model.
4. Residual/Jacobian and diagnostic configuration hashes must explicitly show
   the coupling mode.
5. Do not implement the mode by solving with one configuration and silently
   re-reading the final state with another unrecorded configuration.

### Acceptance

- Residual/Jacobian difference between `postprocess_only` and avalanche-off:
  `<= 1e-12` relative or the relevant absolute near-zero gate.
- Process fields from `postprocess_only` are nonzero and reproduce a separate
  fixed-state evaluation to `<= 1e-12`.
- Existing default-path snapshots are unchanged.

### Outcome

- `vela_postprocess_branch_verified`; or
- `postprocess_residual_leakage` and stop.

Executed outcome: `vela_postprocess_branch_verified`. The opt-in
`postprocess_only` mode is excluded from both coupled and Gummel solver source
assembly while retaining fixed-state avalanche process diagnostics. Focused
tests verify exact residual/Jacobian equality with `model:none`, nonzero
diagnostic source equality with `self_consistent`, configuration validation,
and explicit `coupling_mode` provenance in the release BV audit.

## 8. Work package 4 - create shared solver-used process records

### Problem

The current VTK path contains useful node-reconstructed fields, while source
assembly uses triangle/edge-local quantities. Equal field names therefore do
not always imply equal support or equal numerical operator.

### Files

Create:

- `include/vela/equation/BVProcessProbe.h`;
- `src/equation/BVProcessProbe.cpp`.

Modify:

- `include/vela/equation/AssemblerUtils.h`;
- `include/vela/equation/FixedStateOperatorAudit.h`;
- `src/equation/FixedStateOperatorAudit.cpp`;
- `src/equation/CoupledDDAssembler.cpp`;
- `src/solver/GummelSolver.cpp`;
- `src/simulation/DCSweep.cpp`;
- `CMakeLists.txt`.

Tests:

- `tests/test_fixed_state_operator_audit.cpp`;
- `tests/test_cell_reconstructed_avalanche.cpp`;
- `tests/test_element_edge_gss_laux_avalanche.cpp`;
- add `tests/test_bv_process_probe.cpp` if separation is clearer.

### Record content

Emit solver-used records for each active cell/local edge/local vertex:

- state endpoints and midpoint density;
- E and electron/hole QFP-gradient vectors;
- low-field mobility, high-field drive, final mobility, and limiter;
- directed SG flux and reconstructed/native-labelled current vector;
- alpha_n/alpha_p;
- Gava_n/Gava_p/Gava_total;
- source weight/measure and qG contribution;
- residual row contribution and scatter targets;
- every floor, clamp, sign, truncation, interpolation, or short-circuit branch;
- active-branch fingerprint.

One shared formula owner must supply Real, AD, residual, Jacobian, and
diagnostic paths. The process probe must not reimplement production formulas
in a CSV writer.

### Output semantics cleanup

1. Add explicitly named `NodeReconstructed*` aliases for VTK node-recovered
   current and mobility fields. Preserve existing names temporarily for
   compatibility and document their provenance.
2. Add solver-used cell/edge outputs where the production operator is
   cell/edge based.
3. Deprecate `ElectronIonIntegral`, `HoleIonIntegral`, and `MeanIonIntegral`
   as Sentaurus-equivalent names. Their current local alpha-length
   accumulation should be exposed as `LocalElectronAlphaLengthProxy`, etc.
   Do not compare these proxies with Sentaurus path ionization integrals.

### Acceptance

- The sum of emitted qG contributions matches the assembled source to
  `<= 1e-12`.
- The emitted residual contribution sum matches the assembled carrier blocks
  to `<= 1e-12`.
- Real, AD, and output records share the same configuration and active-branch
  hashes.
- Acute, right, obtuse, reversed-orientation, contact-adjacent, and interior
  tests pass for both carriers.

### Outcome

- `solver_used_process_records_verified`; or
- `diagnostic_formula_drift` and stop.

Executed outcome: `solver_used_process_records_verified`. The normalized
`BVProcessProbe` consumes the existing production SG-edge, triangle-GSS, and
element-edge/GSS-Laux records, emits carrier-split cell/edge/vertex state,
mobility, current, alpha, source, qG, residual-scatter, branch, and
fingerprint data, and closes its source and residual sums to the production
operator at `<= 1e-12`. The coupled assembler exposes the same complete
configuration fingerprint used by the output records. Acute, right, obtuse,
reversed-orientation, contact/interior, both-carrier, and postprocess-only
tests pass. Explicit `NodeReconstructed*` and `Local*AlphaLengthProxy` VTK
aliases document the old fields without removing compatibility names.

## 9. Work package 5 - record every continuation and Newton attempt

### Problem

The current `newton_history` CSV is written only after a sweep point has
converged. This loses the rejected transition that may select the wrong BV
branch.

### Files

Modify:

- `include/vela/solver/NewtonSolver.h`;
- `src/solver/NewtonSolver.cpp`;
- `include/vela/simulation/DCSweep.h`;
- `src/simulation/DCSweep.cpp`;
- `docs/config_schema.md`;
- `configs/schema/vela-simulation.schema.json`.

Tests:

- `tests/test_newton_solver.cpp`;
- `tests/test_dc_sweep.cpp`.

### Required records

`newton_attempts.csv`:

- run/segment/attempt identifiers;
- parent accepted bias and state hash;
- requested and actual target bias;
- initial-state/predictor hash;
- solver and handoff stage;
- accepted/rejected status and typed reason;
- retry number, step size, damping/clamp summary;
- initial/final residual and state hash.

`newton_iterations.csv`:

- attempt and iteration identifiers;
- block residuals and row-scaled residuals;
- raw/applied update norms;
- damping and line-search attempts;
- carrier-row convergence;
- top residual node per equation;
- source/Jacobian active-branch fingerprint.

Both files must include failed and rejected attempts. Existing summary columns
may remain for compatibility.

Also document and expose the already implemented
`write_state_every_point_prefix` in the public config schema and PN2D template
workflow.

### Acceptance

- A synthetic failed transition writes a complete attempt and iteration
  history before retry/shrink.
- The sealed PN2D transition
  `-19.692187499999644 -> -19.693749999999643 V` reproduces as a mandatory
  regression fixture with `max_iterations`, unless a separately authorized
  candidate intentionally changes it and passes the full parity contract.
- A successful retry links to the same parent state and the rejected attempt.
- Repeated runs produce identical typed transition sequences and state hashes.
- Diagnostics do not change convergence or final states.

### Outcome

- `complete_nonlinear_trace_available`; or
- `rejected_transition_unobserved` and stop.

Executed outcome: `complete_nonlinear_trace_available`. The compatibility
accepted-iteration CSV remains unchanged, while the new attempt and iteration
CSVs record rejected transitions before shrink/retry, deterministic parent and
state hashes, requested versus actual targets, typed reasons, raw and
row-scaled block residuals, residual peaks, line-search decisions, and
state-dependent solver-used source/Jacobian active-branch fingerprints. A
synthetic large-step failure proves same-parent retry linkage and byte-stable
repeated traces; a diagnostics-off control produces the same accepted states.
The sealed repository fixture reproduces
`-19.692187499999644 -> -19.693749999999643 V` as `max_iterations` after 40
Newton iterations and emits the initial row plus all 40 iteration rows.

## 10. Work package 6 - close the non-smooth source Jacobian blocker

This work package implements Task 4 of the companion parity plan after work
packages 1, 4, and 5 make the active branch observable.

### Actions

1. Reproduce the accepted `-20 V` source-only probe from an imported
   Sentaurus state or another explicitly fixed assembled state. This probe
   does not require the current self-consistent Vela avalanche-on sweep to
   converge at `-20 V`.
2. Record the exact active branch for QFP differences, mobility drive,
   current magnitude/sign, alpha, source partition, and scatter.
3. Compare shared Real, local forward AD, branch-frozen numerical derivative,
   and symmetric FD over at least three same-branch steps.
4. Determine whether the mismatch is a missing derivative, formula drift,
   cancellation, or a true nondifferentiable branch.
5. For a real branch crossing, use a semismooth/active-set-consistent
   derivative or a branch-separated test. Do not add arbitrary smoothing.

### Acceptance

- Nonzero same-branch maximum analytic/reference difference `<= 1e-8`
  relative.
- Near-zero maximum `<= 1e-12` absolute.
- Residual/Jacobian/process-record hashes agree.
- Contribution sums close independently to `<= 1e-12`.

### Outcome

- `source_jacobian_dependency_identified_and_closed` authorizes work package 7;
- any other companion-plan Task 4 outcome stops nonlinear candidate work.

## 11. Work package 7 - build the paired process-chain analyzer

### Files

Create:

- `scripts/analyze_pn2d_bv_process_chain.py`;
- `tests/regression/test_pn2d_bv_process_chain.py`.

Reuse rather than duplicate:

- `scripts/analyze_pn2d_high_bias_oracle.py`;
- `scripts/analyze_pn2d_high_bias_same_support.py`;
- `scripts/run_pn2d_sentaurus_fixed_state_sg_audit.py`;
- the curve/knee analyzer introduced by Task 2 of the companion plan.

### Analysis sequence

For every exact branch/bias pair, compare in this fixed order:

1. potential and QFP;
2. density and effective intrinsic-density availability;
3. E and QFP-gradient vectors;
4. low/high-field mobility;
5. SG/current-vector stage;
6. alpha;
7. carrier-split generation;
8. geometric qG contribution;
9. assembled residual/Jacobian and first Newton update; and
10. terminal carrier and total current.

For scalars report signed difference, relative difference, and log error where
defined. For vectors report magnitude ratio, angle, component error, and sign.
For spatial support report:

- first/max-error support;
- hotspot and fixed-hotspot values;
- active-support overlap;
- source-weighted centroid;
- cumulative 10/50/90% support; and
- electron/hole, contact/interior, and topology classes.

The analyzer must perform four comparisons:

1. Sentaurus versus Vela avalanche-off;
2. Sentaurus versus Vela IIC/postprocess-only;
3. Sentaurus versus Vela avalanche-on;
4. within each simulator, IIC/postprocess-only versus avalanche-on.

### Outputs

- `stage_summary.csv`;
- `support_summary.csv`;
- `hotspot_chain.csv`;
- `first_departure.json`;
- `source_terminal_closure.csv`;
- `newton_first_update.csv`;
- process-chain and hotspot figures; and
- `acceptance.json` with one typed outcome.

### Acceptance

- Synthetic fixtures recover a deliberately injected error at every stage.
- The earliest stage is invariant to row ordering and unrelated tail values.
- A claimed causal stage is reproduced at two adjacent knee biases.
- `ReadMeasure`/CurrentPlot and Vela assembled-source closure each pass their
  independent gates before cross-simulator conclusions are emitted.

### Outcomes

- `fixed_state_operator_cause`;
- `density_qfp_feedback_cause`;
- `mobility_current_feedback_cause`;
- `source_support_feedback_cause`;
- `contact_boundary_cause`;
- `continuation_solver_path_cause`;
- `proprietary_operator_difference`; or
- `insufficient_observation`.

Only an outcome with two-bias causal evidence may authorize a physics/operator
candidate.

## 12. Work package 8 - implement and qualify one minimal candidate

Follow Tasks 7-11 of the companion parity plan. Each experiment changes one
axis only and keeps defaults unchanged until final review.

Authorized candidate categories are limited to evidence produced by work
package 7:

- a source-Jacobian correction;
- a solver-used source/current support correction;
- a missing high-field mobility feedback derivative;
- a contact-boundary correction; or
- a continuation-only control proving branch invariance.

The candidate must improve both the named internal process metric and the
external curve/knee metrics. A visually better curve without the predicted
internal correction is rejected.

Qualification order:

1. focused unit and derivative tests;
2. fixed-state and first-update tests;
3. duplicate M0 global/knee on/off sweeps;
4. curve, gain, slope, knee, closure, and determinism gates;
5. dose-preserving paired mesh refinement;
6. full Release CTest;
7. 201-point forward-IV regression; and
8. independent scientific and code review.

The candidate remains opt-in unless every companion-plan gate passes.

## 13. Dependency order and parallel-safe work

Mandatory dependency chain:

```text
WP0 -> WP1 -> WP2 -> WP2.1 -> scientific Task 2 verifier
                                -> WP3 -> WP4 -> WP5 -> WP6 -> WP7 -> WP8
```

WP0-WP4 and scientific Task 2 are complete. WP5 is the next gate. After WP1
is stable, these
implementation streams may proceed independently, but their acceptance gates
remain ordered:

- Sentaurus deck/runner implementation in WP2;
- Vela postprocess-only implementation in WP3; and
- initial process-record types/tests in WP4.

Do not run the final paired analyzer until both simulators validate against the
same schema. Do not begin a nonlinear correction before WP6 closes.

## 14. Verification commands

Initialize the required Windows UCRT64 environment:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\vela-tcad"
```

Focused Python tests:

```powershell
python -m unittest `
  tests.regression.test_pn2d_bv_process_contract `
  tests.regression.test_pn2d_bv_process_matrix_vm `
  tests.regression.test_pn2d_bv_process_chain -v
```

Focused C++ tests after the corresponding targets are built:

```powershell
build-release\test_impact_ionization.exe
build-release\test_element_edge_gss_laux_avalanche.exe
build-release\test_cell_reconstructed_avalanche.exe
build-release\test_fixed_state_operator_audit.exe
build-release\test_newton_solver.exe
build-release\test_dc_sweep.exe
```

Full verification:

```powershell
cmake --preset windows-ucrt64-release
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
D:\msys64\usr\bin\git.exe diff --check
```

Live Sentaurus commands must use the existing `sentaurus` SSH alias and keep
all fetched artifacts under run-specific ignored `build-release/` roots.

## 15. Commit boundaries

Use small reviewable commits:

1. common process schema and fail-closed contract tests;
2. Sentaurus IIC/exact-snapshot runner and dry-run tests;
3. Vela postprocess-only coupling mode and focused tests;
4. shared solver-used process records and output-semantics tests;
5. complete Newton/continuation attempt trace;
6. Jacobian blocker correction, only if authorized;
7. paired process-chain analyzer and synthetic tests;
8. minimal opt-in candidate and qualification summaries.

Before every commit:

- inspect `git status --short`;
- stage explicit source/test/document paths only;
- do not stage `tmp/`;
- do not stage `build*/`, TDR, PLT, VTK, CSV, state, log, or generated figure
  artifacts; and
- run `git diff --cached --check`.

## 16. Work-package report template

Every work package must report before the next begins:

1. work-package number and typed outcome;
2. commit and exact commands;
3. input/config/executable/output hashes;
4. requested and actual bias coverage;
5. observed metrics against thresholds;
6. first failed bias/support/carrier/quantity/branch, if any;
7. artifact paths;
8. scoped `git status --short`; and
9. explicit authorization or prohibition for the next work package.

## 17. Next execution slice

WP0-WP4 and scientific Task 2 are complete. The next slice is:

1. implement WP5 without changing a production physics formula; and
2. stop at WP6's authorization gate before any nonlinear correction.

The deterministic current-branch failure is an expected observed condition,
not permission to relax convergence or substitute missing curve rows.
