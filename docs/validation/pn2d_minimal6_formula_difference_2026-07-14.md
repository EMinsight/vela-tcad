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
`build-release/pn2d-minimal6-formula-diff-task5-reviewfix-20260717-a/`.
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
labels, and disclaimer are unclipped. The manifest is marked `reviewed` by
Codex on 2026-07-17. The complete final Task 1-5 Python gate passes `166/166`;
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
