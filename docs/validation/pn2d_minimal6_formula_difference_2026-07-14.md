# PN2D minimal6 formula-difference evidence

## Authoritative Phase A state chain (current)

- Regenerated v1 run: `minimal6_states_v2_20260717_000955`.
- Regenerated v1 manifest SHA-256: `4cc193b9d3ae76ea919d60c6f098902f315e1c0c59de8867b125cbe261769cff`.
- Recovery-validation SHA-256: `f54737dffc06279a53d14cb6a3b82e7d2aa5265875fe29188190ca9b13879d2c`.
- Sealed v2 run: `minimal6_states_v2_sealed_20260717_000955`.
- State identity: `sketch/mirror x 0/-12/-19 V`; all six entries passed and `outputs_complete=true`.
- Input source: approved remote regeneration followed by six hash-bound local UCRT64 operator replays.
- Fixed-state audit: `build-release/reference_tcad/pn2d_sentaurus2018_minimal6/reports/minimal6_fixed_state_audit_task5_20260717_000955/`.

This is the sole current authoritative Phase A chain. The earlier
`b44ad95d...` local-recovery source and `f519aea...` seal are retained only as
historical, superseded evidence. The current immutable sibling seal validates
against `schemas/vela.pn2d_minimal6_states.v2.schema.json`; every state preserves
relative raw-artifact paths and SHA-256 values, and the source-manifest and
raw-member mutation tests fail closed.

The installed Sentaurus Device 2018 deck already requests
`AvalancheGeneration`, `eVelocity`, `hVelocity`, `eIonIntegral`,
`hIonIntegral`, and `MeanIonIntegral`. The validated raw-to-normalized
identities are:

| normalized name | raw name | components | unit | semantic role |
|---|---|---:|---|---|
| `sentaurus_native_avalanche_generation` | `ImpactIonization` | 1 | `cm^-3*s^-1` | `native_avalanche_generation` |
| `sentaurus_electron_speed` | `eVelocity` | 1 | `cm*s^-1` | `carrier_speed` |
| `sentaurus_hole_speed` | `hVelocity` | 1 | `cm*s^-1` | `carrier_speed` |
| `sentaurus_electron_ionization_integral` | `eIonIntegral` | 1 | `1` | `path_ionization_integral` |
| `sentaurus_hole_ionization_integral` | `hIonIntegral` | 1 | `1` | `path_ionization_integral` |
| `sentaurus_mean_ionization_integral` | `MeanIonIntegral` | 1 | `1` | `path_ionization_integral` |

The C++ fixed-state replay validated all six immutable provenance records,
including full argv arrays and executable/config/input/output hashes. The
joined audit passed with `36/54/24` node/edge/triangle rows, maximum imported
state parity `0.0` (gate `<1e-12`), and maximum C++/Python formula error
`1.79954895429546e-12` (gate `<5e-12`).

## Scope

This is a diagnostic evidence package. It does not establish a physical BV curve or alter production solver defaults.
## Task 5 figure and documentation contract

The formula-difference CLI emits the fixed PNG/PDF pairs `gradient`, `current_alpha`, `source_waterfall`, `interaction`, and `topology_symmetry`, plus `figure_manifest.json`. The manifest retains units, the diagnostic disclaimer, and a reviewer/date/checklist entry. A visual inspection remains a Phase A gate; an automatically produced file is not evidence of a physical BV curve.

The reports maintain these source labels without substitution: native Sentaurus `ImpactIonization`, the separately labelled Sentaurus `alpha*|J|/q` reconstruction, and the separately labelled Vela `alpha*flux*partial_volume` reconstruction. Producer-owned pair classifications distinguish explicit geometric zeros from exact available zeros and unavailable values. Geometric zeros are marked without log flooring; unavailable values are shown as `N/A` and are not plotted as numeric zeros.

## Root-cause implementation map

| Classified factor | Sentaurus deck or parameter entry | Independent diagnostic implementation | Exact production C++ entry |
| --- | --- | --- | --- |
| `ni_eff/BGN` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) `EffectiveIntrinsicDensity(OldSlotboom)` and [models.par](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par) `OldSlotboom` | [physics.py](../../scripts/pn2d_minimal6_diagnostics/physics.py) `infer_ni_eff`; [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_ni` | [ScharfetterGummel.h](../../include/vela/discretization/ScharfetterGummel.h) variable-`ni` electron/hole SG helpers |
| `gradient_recovery` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) exported potential/current controls | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_gradient`; [support.py](../../scripts/pn2d_minimal6_diagnostics/support.py) projections | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `cellScalarGradient` / `edgeAveragedCellScalarGradient` |
| `mobility` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) `Mobility(DopingDependence HighFieldSaturation)` and [models.par](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par) mobility entries | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_mobility` | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `edgeMobility` / `triangleGssEndpointAveragedMobility` |
| `current_semantics` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) `eCurrentDensity/Vector` / `hCurrentDensity/Vector` | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_current` | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `selectAvalancheCurrentFluxProxy` / `sgEdgeCurrentAvalancheSourceRecords` |
| `impact_driving_field` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) `Avalanche(VanOverstraeten)` | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_impact_field` | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `electronAvalancheDrivingField` / `holeAvalancheDrivingField` |
| `alpha_law` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) alpha outputs and [models.par](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par) `vanOverstraetendeMan` | [physics.py](../../scripts/pn2d_minimal6_diagnostics/physics.py) `van_overstraeten_alpha`; [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_alpha` | [ImpactIonizationModel.cpp](../../src/physics/ImpactIonizationModel.cpp) `VanOverstraetenImpactIonization::{electron,hole}Coefficient` |
| `partial_volume` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) native `AvalancheGeneration` support control | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_partial_volume`; [support.py](../../scripts/pn2d_minimal6_diagnostics/support.py) `integrate_cell_field` | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `avalancheSourceEdgeArea` / `sgEdgeCurrentAvalancheSourceRecords` |
| `source_to_node_mapping` | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) native `AvalancheGeneration` mapping boundary | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py) `_formula_source_mapping`; [support.py](../../scripts/pn2d_minimal6_diagnostics/support.py) `local_edge_sources_to_nodes` | [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) `sgEdgeCurrentAvalancheSourceComponentIntegrals`; [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp) continuity-source assembly |

Parameter agreement is a control, not a causal conclusion. The current formula report remains `insufficient_data` until raw Vela state/operator inputs permit declared counterfactual substitutions and closure.

## Task 5 re-execution evidence (2026-07-17)

Task 5 was re-executed from clean baseline `f1d8c67`. The current TDD cycle
corrected the `current_alpha` figure to use the ledger's SI units (`A/m^2` and
`m^-1`) and to render signed electron/hole current projections plus distinct
electron/hole alpha series. The generated Markdown now links all eight
classified factors to tracked Sentaurus, independent Python, and exact C++
entries. Manual QA then found the 24-label interaction bar chart unreadable;
a second RED/GREEN cycle replaced it with a 6-state by 4-ordered-pair annotated
heatmap without dropping any interaction record. The failed bar-chart package
is preserved under
`build-release/pn2d-minimal6-formula-diff-task5-real-20260717_000955-a/`.

The first approved regeneration run,
`minimal6_states_v2_20260716_234905`, is preserved with `5/6` passed states and
`outputs_complete=false` after the 15-minute local orchestration timeout. It
was not overwritten or sealed. A distinct unchanged-physics regeneration,
`minimal6_states_v2_20260717_000955`, completed all six exact states. The
importer omitted only the local replay inputs `mesh.json` and `audit.json`;
their restored copies were byte/hash checked against the tracked configuration
and the previously validated real candidate, while all regenerated node and
element hashes matched the prior six-state topology.

Six actual local UCRT64 operator replays were recorded and hash-bound before
the fixed-state audit. The regenerated v1 manifest SHA-256 is
`4cc193b9d3ae76ea919d60c6f098902f315e1c0c59de8867b125cbe261769cff`;
its recovery-validation SHA-256 is
`f54737dffc06279a53d14cb6a3b82e7d2aa5265875fe29188190ca9b13879d2c`.
The immutable sibling seal
`minimal6_states_v2_sealed_20260717_000955` validates as
`vela.pn2d_minimal6_states.v2` with six states. The bound audit under
`build-release/reference_tcad/pn2d_sentaurus2018_minimal6/reports/minimal6_fixed_state_audit_task5_20260717_000955/`
passes exact `36/54/24` identities, state parity `0.0`, C++/Python formula
hybrid error `1.7995489542954602e-12`, and replay provenance.

Two fresh real reports are byte-identical for every deterministic report
artifact. The reviewed package is
`build-release/pn2d-minimal6-formula-diff-task5-symmetryfix-20260717-a/`.
It validates against `vela.pn2d_minimal6_formula_difference.v1`, reports
`insufficient_data` because `ni_eff/BGN`, `impact_driving_field`, `alpha_law`,
and `source_to_node_mapping` remain unavailable, and emits no
`dominant_factor`. Both forward/reverse closure error and named-residual
mismatch are `0.0 dex`; all 24 interaction records are retained.

| artifact | SHA-256 |
|---|---|
| `root_cause_summary.json` | `3229cf3c504215db88e1600ea50636bdf363d96e7aad553e871eac031c292493` |
| `quantity_ledger.csv` | `38aa4dc6623f944190be5fd575e7a87da2e8fa7532a9dd45288bbe809dca671f` |
| `factor_waterfall.csv` | `2df94955bee1c11394fbaad7be67742cbe679468e37538941ad7f24cb652aca0` |
| `root_cause_summary.md` | `86bb168578a6467747b87f88cf86534a6db559b9d0b21c5b8dc2ff1ab3d0c335` |
| `figure_manifest.json` | `4776c0a91cf42b53213e835257c4142e14fd595bc5f5ac14ce47b4eefcf08d35` |

Manual QA passed all five PNG/PDF pairs: exact state identities, sketch/mirror
symmetry, signed gradients and carrier currents, SI units, explicit zero
handling, three distinct native/reconstructed source labels, waterfall
closure, all 24 readable heatmap cells, and the mandatory disclaimer. All 18
real source rows are explicitly `available`; independent visual QA confirmed
that no false geometric-zero or `N/A` marker is present and that the color key,
labels, and disclaimer are unclipped. The topology-symmetry PNG remains
byte-identical to the previously reviewed real figure (SHA-256
`15c1c3b9a6b9f4166215094e5f0839edd23636ee8103cec49c9598d32f6572c3`);
the real all-available matrix correctly requires no gap markers. The manifest
is marked `reviewed` by
Codex on 2026-07-17. The complete final Task 1-5 Python gate passes `168/168`;
the release build exits `0`; C++ controls
pass `511` assertions / `40` cases and `83` assertions / `14` cases.

Phase A is complete with a typed, integrity-valid `insufficient_data` result.
Under the plan's terminal-state rules, Phase B may proceed without treating
this result as a dominant-factor claim or a physical BV conclusion.

## Task 8 final evidence (2026-07-15)

Answer first: the regenerated real six-state report validates, but the strict
root-cause result remains `insufficient_data`. Raw Sentaurus ledgers and the
production C++ Vela replay ledgers are complete; the raw Vela state required
for the declared named counterfactual substitutions is not available. No
factor ranking or residual closure was fabricated.

The reviewed report is under
`build-release/pn2d-minimal6-formula-diff-task8-20260715/`. Its schema is
`vela.pn2d_minimal6_formula_difference.v1`, and it contains exact
`36 / 54 / 24` node/edge/triangle identities. The five required figure pairs
(`gradient`, `current_alpha`, `source_waterfall`, `interaction`, and
`topology_symmetry`) passed manual label, unit, disclaimer, zero/unavailable,
and clipping checks. The figure manifest is marked `reviewed` by Codex on
2026-07-15.

| artifact | SHA-256 |
|---|---|
| `root_cause_summary.json` | `689ae18c4c608c850bb10043c9bebc7bf68f150481aa80a0cfabce91e1f7f158` |
| `quantity_ledger.csv` | `2a0a36dc22f4589d33fd2133718e40c6c049f2fd1e9452605cdce09fdfbda9d9` |
| `figure_manifest.json` | `7a8ad92345b94e6ba1b1adc2963b779f4e547a53287f1d4f78360b91800d7700` |

The independent nonlinear evidence is retained separately:

- Vela: `build-release/pn2d-minimal6-vela-task8-final-r3-20260715/`;
  sketch and mirror accept exact `-1 V`, close their production SG native and
  reconstructed source totals exactly, then both reject `-1 -> -2 V` with
  exit code 1 and `nonfinite_residual`.
- Sentaurus:
  `build-release/pn2d-minimal6-sentaurus-task8-final-r2-20260715/`; both
  topologies accept every exact checkpoint from `0` through `-20 V` (42 rows),
  with no rejected transition. The clean `0 V` rows pass the scalar/vector
  field-contract selector.
- Comparison:
  `build-release/pn2d-minimal6-comparison-task8-final-r2-20260715/`; the only
  common exact bias is `-1 V`. Vela/Sentaurus terminal-current ratios are
  `1.4457180849345512e-4` (sketch) and `1.5145723823453054e-4` (mirror).
  Under the recorded `v1` threshold (`ratio <= 1e-3`), both are
  leakage-like. Counterfactual gap closure remains unidentifiable.

All of these curves retain the mandatory statement:

> minimal6 diagnostic sweep; not a physical BV curve

## Task 8 fresh evidence and nonlinear recheck (2026-07-17)

The fresh Task 8 fixed-state package is
`build-release/pn2d-minimal6-formula-diff-task8-20260717-a/`. It is
schema-valid and deterministic: `root_cause_summary.json`,
`quantity_ledger.csv`, `factor_waterfall.csv`, `root_cause_summary.md`, and
`figure_manifest.json` are byte-identical to the final Task 5 package. The
report SHA-256 is
`3229cf3c504215db88e1600ea50636bdf363d96e7aad553e871eac031c292493`;
the figure-manifest SHA-256 is
`4776c0a91cf42b53213e835257c4142e14fd595bc5f5ac14ce47b4eefcf08d35`.
It retains six states, 24 interactions, no dominant factor, and typed
`insufficient_data`; both forward and reverse waterfall closure mismatch are
`0.0 dex`.

The fresh nonlinear comparison is
`build-release/pn2d-minimal6-comparison-task8-fresh-20260717-a/`, JSON SHA-256
`2d8b2d86f119964950ff6a620e252a9509bac9edd82cd0140a5df1bcb2d56ea0`.
It binds 40 exact common `-1..-20 V` checkpoints. The fixed-state recheck at
`0 V` is `unidentifiable` because no common Vela checkpoint exists. At
`-12 V` and `-19 V`, four hash-addressed nonlinear states exist per bias, but
the verified nonlinear ledger-input bundle required to rerun the Task 4
counterfactual chain is absent; both rankings therefore remain typed
`unidentifiable`. No embedded checkpoint summary is trusted as a substitute.

All 156 eligible sweep log gaps are retained as observations with
`decomposition_status=unidentifiable`, empty named contributions, and null
closure error; no tautological decomposition is claimed. All 40 branch rows
are `unidentified`, and no multiplication branch or physical BV is inferred.
The full hash audit SHA-256 is
`73650d0625fb7a49af6581c70b48afcd6e50b0536857ea8bad738dd441cf512f`.

> minimal6 diagnostic sweep; not a physical BV curve

## Task 9 physics inverse audit (2026-07-22)

Answer first: the exact 40-state inverse audit is deterministic and
integrity-valid. It identifies six typed candidate operators, rejects two
direct quasi-Fermi-gradient generation candidates, and classifies seven
current-related candidates as `confounded`. It does **not** identify a complete
replacement formula: the seven-stage replacement chain stops with
`status=missing_field` at the unavailable `mobility` factor, so its baseline,
full replacement, and closure are all `null`. No production C++ formula was
changed.

The authoritative report is
`build-release/pn2d-minimal6-physics-inverse-audit-20260722-g/`; the independent
reproduction is the sibling root ending in `-h`. Both contain 25 files with
identical relative paths and identical SHA-256 values for all 25 files. Each
independent verifier passes all 15 checks over 22 report artifacts and 1,274
input members. The report-manifest SHA-256 is
`0e979e1ae323d0413a16f4f87879271285e188ad1a523a2f201ccb1b21f3193c`,
the independently recomputed scientific-payload SHA-256 is
`e1e4edddca233958145ecfc9b70d6a14918372b85a77bdd832789d8d81ffed6b`,
and the canonical input-manifest SHA-256 is
`1bf13239725fbdb9439b7d20164a0611cd601076cbc8c1ec3520da522fa71300`.
The phase base is `a5524cf`, `production_cpp_changed=false`, and the unique
Sentaurus version is `O-2018.06-SP2`.

The earlier `-c/-d` roots are retained only as historical pre-effective-deck
packages. The `-e/-f` roots are diagnostic failed-generation roots that exposed
an exact Python binary64 operation-order mismatch in the independent verifier;
although `-e` passes after that verifier repair, neither root is promoted to
authority. Fresh `-g/-h` generation and post-generation independent verification
both pass without tolerance-based comparison.

The report and independent verifier also enforce the labelled vertical mirror
contract before accepting any inverse evidence: coordinates transform as
`(x,y)->(x,0.5 um-y)`, vector components as `(vx,vy)->(vx,-vy)`, and the
zero-based node map is `0<->4`, `1<->5`, `2<->3`. With relative tolerance
`1e-9` and a dimensionally scoped `1e-8 V/m` absolute tolerance, the sealed
inputs contain 3,240 valid mirror pairs and 1,080 matching nonvalid pairs, with
zero mismatches and zero unpaired samples. The earlier `-a/-b` roots are
retained only as provisional, pre-mirror-gate packages and are not
authoritative.

The `production_cpp_changed=false` payload flag is evidence metadata rather
than an internally executed Git check. The required external phase-base diff
against `a5524cf` is therefore retained as the fail-closed production-source
guard; it confirms no changes under `include/` or `src/`.

All 40 hash-bound Vela decks resolve to the same complete 38-key canonical
effective impact-ionization configuration. The effective values are
`model=van_overstraeten`, `parameter_set=default`,
`driving_force=quasi_fermi_gradient`, `generation=current_density`,
`current_approximation=density_gradient`, `A_scale=1.0`, and `B_scale=1.0`;
the derived thermal voltage is `0.025851999786435 V`. Both the report builder
and independent verifier parse and validate this effective deck configuration
rather than assuming header defaults.

The remote supplemental root completed exactly 40/40 requested sketch/mirror
states from `-1` through `-20 V`, with zero prepared or failed states and
`outputs_complete=true`. Its manifest SHA-256 is
`675f3f42fbe463700330c630ac056bd5270cb240b7f9b4ca1378412dc2e3433d`.
The 1,376 members belonging to the 32 states that predated resume remained
byte-identical before and after resume; their aggregate SHA-256 stayed
`5d972825781b4bf8273390f88f5f05b24cce74bc3d7f6549e04e9587aa323c15`.

The three canonical input roots are independently sealed:

| root | manifest SHA-256 | seal SHA-256 |
|---|---|---|
| `vela` | `131b4fe53fefa90d6c43ede6317c9a24a78111ecaa1d94dc664a4f4af3d4b883` | `0b0edea45f07d39adba38f232bab02e44ca423066012032aebba1c80307a4d5c` |
| `sentaurus` | `7ed91f79720211ab943eb22e2a8939645cea940651a15cbc50ec0ca05ce0cfb2` | `94f82c24ca137965bbffa4756908ad9750f31bf104e74e4902e5e2aab8364041` |
| `supplemental` | `eeb72a54c48b4e56258075f50bbdbddf3a1545ba582537e79f63ab767b17eb59` | `7b2679bb840edbd9f772f78f46b2fd408f05af2cbb058d9c41ca5e1175822cff` |

The report contains 8,640 canonical observations over seven discovery and 33
holdout states: 6,480 are `valid` and 2,160 are explicitly
`missing_field`. Vela contributes all 2,160 gaps: 480 ElectricField
components, 480 electron-current components, 480 hole-current components,
240 ImpactIonization values, 240 electron-alpha values, and 240 hole-alpha
values. These gaps are never replaced with projected edge/cell values or
zeros.

The report contains 90 metric rows and 15 final candidate classifications:

- `identified`: `node_area_weighted_minus_grad_psi`,
  `edge_area_weighted_minus_grad_psi`, `triangle_minus_grad_psi`,
  `electric_field_magnitude`, `electric_field_current_aligned`, and
  `signed_edge_minus_delta_psi_over_h`;
- `rejected`: `qf_gradient_magnitude` and
  `qf_gradient_current_aligned`;
- `confounded`: `signed_edge_sg_density_current`,
  `signed_edge_drift_diffusion_current`,
  `node_area_weighted_qf_gradient_current`,
  `current_inverted_qf_gradient`, `signed_edge_qf_difference_current`,
  `edge_area_weighted_qf_gradient_current`, and
  `triangle_qf_gradient_current`.

`identified` here means that the declared candidate passes discovery, holdout,
and combined numerical gates; it is not a claim that the complete production
operator has been recovered. For the combined split, the three `-grad(psi)`
reconstructions have zero median direction error and median relative magnitude
errors from `1.8807661440764189e-16` to `4.6797221468914373e-13`.
The electric-field generation candidates have combined integrated-generation
median absolute errors of `5.251003410267687e-05 dex` (magnitude) and
`5.2286790636912786e-05 dex` (current-aligned). In contrast, the two direct
quasi-Fermi-gradient generation candidates have combined integrated-generation
median errors near `0.745609 dex` and fail the acceptance gate.

The current-inversion family remains confounded in discovery, holdout, and
combined splits because mobility is not independently available for the full
cross-solver substitution. For example, `current_inverted_qf_gradient` has a
combined transport median error of `0.21694243743356303 dex` and median
direction error of `89.741532388140698 degrees`, but these values cannot isolate
current semantics from mobility and quasi-Fermi-gradient recovery.

The seven-stage replacement matrix is therefore a deterministic typed QA
control, not a causal substitution into production physics. Its dependency
order is gradient recovery, mobility, current semantics, impact driving field,
alpha law, geometric integration, and source-to-node mapping. Because the
second factor is unavailable, there are no one-factor, forward, reverse, or
adjacent-interaction sequences, and no numerical closure is claimed. The next
formula-identification experiment requires native same-support Vela electric
field and signed current vectors, independently comparable carrier mobility
and quasi-Fermi gradients, alpha values, and volumetric plus integrated
generation. No production formula replacement is justified by the present
audit.

> minimal6 inverse audit; diagnostic identities and arithmetic closure are not a physical BV curve or a production-formula change
