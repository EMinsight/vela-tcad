# PN2D Minimal6 current-alignment Phase B-G revised plan

Date: 2026-07-24

Status: Phase F complete with typed `model_difference`; Phase G is next.

## Objective

Determine why the self-consistent Vela electron and hole quasi-Fermi
potentials differ from Sentaurus, then reduce the resulting carrier-density,
mobility, current, and avalanche-source differences without changing a
production formula that has not been causally falsified.

This plan supersedes any Phase B-G action that treats the SG current algebra,
carrier sign, or Minimal6 box geometry as the leading suspected defect.

## Authoritative evidence baseline

The following results are frozen inputs to this plan:

- Vela and Sentaurus internal electrostatic potential differ by at most
  `2.24802e-11 V` at nodes 1 and 5.
- Independently solved internal electron/hole QFP differs by approximately
  `0.33 V`, with maxima near `0.368 V`.
- With Sentaurus electrostatic and QFP potentials imported, Vela BGN
  recomputes all 240 node carrier densities within `4.426181e-6 dex`.
- The earlier Vela-gradient element mobility control differs by medians
  `0.052688 dex` for electrons and `0.047814 dex` for holes. Phase D showed
  that this is not a native-Sentaurus-field comparison: using the exported
  native element QFP gradient gives medians `0.593706 dex` and `0.047364 dex`.
- On the box-edge support after imported state and recomputed Vela mobility,
  median mobility/current differences are approximately `0.040736 dex`
  electron and `0.039068 dex` hole. The high-bias central-edge tail reaches
  `1.063346 dex` and `0.710159 dex`, respectively, on currents near
  `1e-20 A/um`.
- The reconstructed Sentaurus box current closes total terminal current to
  `1.137692e-7` relative and internal total-current KCL to `3.163255e-9`.
- Vela production SG current is replayed by the same box algebra to
  `5.295713e-14` maximum relative difference.
- Vela and Sentaurus box/cotangent geometry is identical on the observed
  Minimal6 support.
- The `sketch` Sentaurus region-cell mapping is
  `0->0, 1->3, 2->2, 3->1` in Vela triangle order. Treating region-cell order
  as triangle-id order is forbidden.

## Revised causal model

The working dependency chain is:

`QFP nonlinear state`
`-> carrier density`
`-> mobility evaluated on a defined support`
`-> SG/box edge current`
`-> continuity residual and avalanche source`
`-> self-consistent branch and terminal current`.

The current evidence supports the following ranking:

1. internal QFP state/continuity equation;
2. mobility parameter and element-to-edge support semantics;
3. source terms, boundary conditions, and nonlinear branch recovery;
4. SG algebra and geometry only as regression guards, not active suspects.

## Phase matrix

| Phase | Revised purpose | Status | Required exit |
|---|---|---|---|
| B | Establish a terminal-closed Sentaurus box-current reference. | Complete | Exact topology, coefficient, terminal-current, KCL, sign, and independent-verification gates pass. |
| C | Run QFP -> density -> mobility -> geometry replacement on the exact 40-state lattice. | Complete | All 40 states and 400 nonzero carrier-edge values are typed; density and full-replacement closures pass; A/B roots are deterministic. |
| D | Align mobility model parameters and element/edge support without fitting current. | Complete | Typed `proprietary_model_difference` plus `support_mismatch`; documented parameter mismatch rejected. |
| E | Localize the QFP continuity residual and nonlinear branch drift. | Complete, including post-fix rerun | Shared `1e-8` source conversion validated; maximum residual `2.107970e-8`; Jacobian relative difference `2.989497e-9`; full physics converges at `-20 V`. |
| F | Re-run self-consistent 40-state sweeps and compare the complete physical chain. | Complete, `model_difference` | 40 exact states converge; electrostatic potential passes; electron QFP is the first failed metric; A/B evidence is deterministic. |
| G | Decide whether a production change is justified; validate and review it if so. | Source-unit patch implemented; final review after F | Retain the minimal source-unit patch and complete full validation and review after the Phase F comparison. |

## Phase B - Sentaurus box-current reference

### Completed work

1. Probed native directed-edge output and recorded
   `schema-valid insufficient_data`.
2. Enumerated 9 global edges, 4 elements, 12 element-local edges, and
   `ReadCoefficient`.
3. Replayed the documented electron QFP-plus and hole QFP-minus box branches.
4. Verified contact electron/hole currents, total-current KCL, and geometric
   zero edges.
5. Rejected element-current-vector projection and wrong QFP-sign controls.

### Frozen gates

- The reference must always be labeled `Sentaurus operator reconstruction`,
  not native directed-edge observation.
- Total terminal current relative error must remain below `2e-7`.
- Internal total-current KCL relative error must remain below `1e-8`.
- Zero coefficient edges must remain exact typed zeros.
- Any new topology must obtain a field-verified region-cell mapping before
  element values are used.

No further native-edge export work is required unless a new vendor API or
documented output location becomes available.

## Phase C - 40-state staged replacement

### Completed work

1. Replayed `mirror/sketch x -1..-20 V`.
2. Corrected `sketch` region-cell order using native element electric field.
3. Replaced QFP, density, element mobility, and geometry in strict order.
4. Retained a Vela-BGN recomputed-density control.
5. Produced two independent byte-identical result roots.

### Frozen interpretation

- QFP replacement fixes current direction but not magnitude.
- Density consistent with imported potentials removes about
  `5.24-5.34 dex` and is the dominant replacement step.
- Mobility removes the remaining approximately `0.04-0.06 dex` median.
- Geometry contributes exactly zero.
- A large relative error on the high-bias central edge must be reported
  together with its approximately `1e-20 A/um` absolute current.

Phase C must be rerun as a regression after any Phase D-G code or parameter
change.

## Phase D - mobility parameter and support alignment

### D1. Freeze the mobility configuration

- Extract the active Sentaurus mobility models and numeric parameters from
  the command and parameter files.
- Record Masetti parameters, high-field model, saturation velocities,
  temperature, carrier statistics, and any interpolation settings.
- Record the corresponding Vela configuration after unit conversion.
- Fail closed on an implicit default or an unsealed parameter source.

### D2. Compare on native element support

- Evaluate the Vela low-field and high-field mobility on each of the 160
  native Sentaurus elements.
- Use the same element QFP-gradient magnitude and explicit doping controls:
  cell-average net doping, arithmetic node mobility, and any documented
  Sentaurus interpolation.
- Preserve low-field, high-field ratio, and final mobility as separate
  columns.
- Retain the already passing velocity-unit parity branch as a regression.

### D3. Construct an explicit box-edge mobility

- For each global edge compute the Sentaurus-equivalent coefficient-weighted
  mobility
  `sum(kappa_element_edge * mobility_element) / sum(kappa_element_edge)`.
- Compare it with the Vela production global-edge mobility and the three
  cell-local Vela mobilities.
- Keep zero-coefficient edges typed and exclude them from dex statistics.
- Do not average element indices until the field-derived cell permutation has
  passed.

### D4. Localize the central-edge tail

- Decompose the 1-5 edge by bias, topology, carrier, adjacent element,
  low-field mobility, QFP-gradient field, saturation factor, and current
  magnitude.
- Determine whether the high-bias tail is caused by one global Vela edge
  mobility replacing two Sentaurus element mobilities, by a parameter
  difference, or by both.
- Include absolute-current weighting in addition to unweighted dex.

### D5. Parameter-substitution controls

- Run no-fit substitutions using only documented Sentaurus parameters.
- A fitted parameter is allowed only as a diagnostic and must not be proposed
  for production.
- Separate `parameter_mismatch`, `support_mismatch`, and
  `proprietary_model_difference`.

### Phase D exit gate

Phase D passes when all 320 native carrier-element samples and all active
box-edge samples have explicit common support and one of these typed outcomes:

1. a documented parameter/configuration mismatch;
2. a documented element-to-edge support mismatch; or
3. an irreducible model difference with bounded residuals.

Targets for a parity candidate are median mobility error at or below
`0.03 dex` and p95 at or below `0.10 dex` on native element support. Missing
the target is not permission to fit; it produces a typed model-difference
result.

### Completed outcome

- All 22 documented numeric carrier parameters match after SI conversion;
  the documented-parameter substitution is an exact no-op.
- On the native Sentaurus element QFP-gradient support, electron mobility has
  median/P95 errors `0.593706/0.703886 dex`; hole mobility has
  `0.047364/0.127339 dex`.
- The coefficient-weighted element-to-global-edge signed-log decomposition
  closes to `1.3877787807814457e-16 dex`.
- The central edge `1-5` tail contains both a native element
  model/interpolation residual and a global-edge support residual; its
  approximately `1e-20 A/um` absolute current remains explicitly reported.
- The no-fit inferred-low-field replay closes to
  `3.7683078942700315e-16` maximum relative error but is diagnostic-only.
- Primary typed result: `proprietary_model_difference`.
- Secondary typed result: `support_mismatch`.
- No production mobility formula change is justified.

Evidence:
`docs/validation/pn2d_minimal6_phase_d_mobility_support_alignment_2026-07-24.md`.

## Phase E - QFP continuity residual and branch localization

### E1. Assemble the imported-state residual

At every exact state, import Sentaurus `psi`, electron QFP, and hole QFP, then:

1. recompute `n,p` with Vela BGN;
2. evaluate three mobility branches:
   Vela production, Sentaurus-equivalent box-edge, and a constant-mobility
   control;
3. compute directed electron/hole SG current;
4. compute SRH and avalanche sources; and
5. assemble the electron and hole continuity residual at nodes 1 and 5.

All terms must use the same sign, volume, depth, and node/edge orientation
contract.

### E2. Produce a residual waterfall

For each carrier/node/state report:

- contact/boundary flux;
- SG divergence;
- SRH contribution;
- impact-ionization contribution;
- scaling/normalization;
- final residual;
- residual change from each mobility branch.

Absolute SI/internal residuals and normalized residuals are both required.

### E3. Audit boundary and Jacobian behavior

- Verify contact QFP Dirichlet values and carrier-density boundary conversion.
- Compare analytic Jacobian entries with finite differences at the imported
  state.
- Record the first Newton/Gummel update from the imported state.
- Distinguish a nonzero physical residual from a correct residual with an
  incorrect Jacobian/update.

### E4. Run controlled branch experiments

Use this order:

1. fixed `psi`, solve electron/hole QFP only;
2. coupled `psi/QFP`, avalanche disabled;
3. coupled `psi/QFP`, SRH disabled control;
4. full configured physics;
5. homotopy from imported QFP to the Vela branch.

Do not tune tolerances or clamp QFP until the unmodified residual and update
are recorded.

### E5. Localize the first divergence

- Run from -1 V upward in magnitude with exact checkpoint identity.
- Identify the first bias/node/carrier where the imported-state residual or
  first nonlinear update departs materially.
- Compare mirror/sketch only after applying their verified support mapping.

### Phase E exit gate

Phase E must identify whether the approximately `0.33 V` QFP difference
originates in:

- mobility coefficients in the continuity equation;
- SRH or avalanche source;
- contact/boundary state;
- residual scaling;
- Jacobian/update implementation;
- nonlinear continuation/branch selection; or
- an explicit Sentaurus/Vela model difference.

No production solver change is allowed until the fixed imported-state
residual and the first update independently reproduce the claimed cause.

## Phase F - self-consistent 40-state comparison

### Entry conditions

- Phase D has a typed mobility outcome.
- Phase E has a fixed-state causal result.
- Any candidate change has focused RED/GREEN tests.
- Phase B/C regression gates still pass.

### Execution

1. Run immutable `mirror` and `sketch` sweeps from 0 to -20 V.
2. Preserve every exact integer checkpoint and the first failed nonlinear
   attempt.
3. Classify the terminal-current branch; convergence alone is insufficient.
4. Compare, in dependency order:
   `psi -> QFP -> n/p -> mobility -> edge current -> terminal current ->
   alpha/source`.
5. Generate baseline and candidate roots independently twice.

### Provisional quantitative targets

These targets must be frozen before the real candidate sweep:

| Quantity | Support | Provisional target |
|---|---|---:|
| internal electrostatic potential | nodes 1/5 | maximum <= `1e-6 V` |
| electron/hole QFP | nodes 1/5 | median <= `0.01 V`, p95 <= `0.025 V` |
| electron/hole density | nodes 1/5 | median <= `0.10 dex`, p95 <= `0.25 dex` |
| mobility | matched native element | median <= `0.05 dex`, p95 <= `0.20 dex` |
| directed current | active reconstructed box edges | median <= `0.10 dex`, p95 <= `0.25 dex`, sign = 100% |
| total terminal current | exact states | median <= `0.10 dex`, no branch-order failure |
| integrated avalanche source | exact states | median <= `0.30 dex` |

Thresholds may be revised only by a pre-run written scientific rationale, not
after observing candidate results.

### Phase F exit outcomes

- `parity_passed`;
- `solver_first_failure` with preserved exact common prefix;
- `model_difference` with the first failed metric and support; or
- `insufficient_data` caused by an integrity or availability limitation.

### Completed outcome

- Both independent sweeps accept all 40 Vela checkpoints with no failed
  transition.
- All 40 states are `multiplication_like`; terminal-current branch order is
  preserved.
- Internal electrostatic potential passes with maximum `2.36565e-11 V`.
- Electron QFP is the first failed dependency metric, with median/P95 errors
  `0.0548737/0.0737589 V`.
- Hole QFP median/P95 errors are `0.0620671/0.0821749 V`.
- Density, native-element mobility on the differing self-consistent state,
  directed current, terminal current, and impact source fail downstream.
- The carrier-edge sum reproduces the Vela Anode current within
  `3.66527e-16` relative, excluding a current integration/sign defect.
- Primary typed result: `model_difference`.
- No additional production formula change is justified.

Evidence:
`docs/validation/pn2d_minimal6_phase_f_self_consistent_2026-07-24.md`.

## Phase G - production decision, regression, and review

### Decision matrix

| Phase D/E result | Allowed action |
|---|---|
| Unit or implementation defect isolated in Vela | Minimal production patch with RED/GREEN test. |
| Sentaurus parameter differs from Vela configuration | Change/configure the comparison deck; do not rewrite the formula. |
| Element/edge support differs | Keep a diagnostic reconstruction or explicitly redesign the discretization; do not fit mobility. |
| QFP residual/Jacobian defect isolated | Patch only the affected residual, boundary, scaling, or derivative. |
| Nonlinear branch/continuation difference | Change continuation only after exact fixed-state residuals agree. |
| Proprietary or undocumented model difference | Record `model_difference`; no production patch. |

### Required validation

- Focused mobility, SG, residual/Jacobian, boundary, and sweep tests.
- Phase B contact/KCL replay.
- Phase C 40-state replacement replay.
- Two byte-identical full candidate roots.
- Full Release build and CTest.
- `ascii_sources`, schema, hash, and diff checks.
- Independent scientific review of units, signs, supports, and conclusions.
- Independent code review of every production change.

### Final deliverables

- revised validation report with a dependency-chain table;
- raw and summarized fixed-state residual evidence;
- full self-consistent sweep evidence or first-failure record;
- formula-decision ledger stating why each changed or unchanged formula is
  justified;
- hash-addressed manifests and independent verification;
- scoped commit(s) that preserve unrelated working-tree changes.

## Global stop conditions

Stop the downstream phase and preserve evidence when any condition occurs:

1. state identity, topology, bias, units, or source hashes do not match;
2. region-cell mapping is missing or fails the electric-field gate;
3. Vela production current cannot be replayed within `1e-12` relative;
4. Sentaurus reconstructed total terminal current exceeds `2e-7` relative or
   KCL exceeds `1e-8`;
5. imported-state density differs by more than `1e-4 dex` before mobility is
   evaluated;
6. a zero/reference-missing value is converted into an invented dex value;
7. a production formula is changed before a fixed-state causal defect is
   independently reproduced;
8. mirror/sketch differ after support mapping without a localized
   topology-specific cause; or
9. a nonlinear run changes configuration, silently skips a checkpoint, or
   loses the first failure.

## Immediate next action

Begin Phase G. Retain the minimal source-unit production patch, run the full
Release/CTest and ASCII/schema/hash/diff checks, prepare the formula-decision
ledger, and request independent scientific and code review. Do not change the
mobility, SG, Poisson, or impact formulas without a new fixed-state causal
falsification.

Phase F evidence: docs/validation/pn2d_minimal6_phase_f_self_consistent_2026-07-24.md.
