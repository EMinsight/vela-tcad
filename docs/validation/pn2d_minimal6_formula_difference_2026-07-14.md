# PN2D minimal6 formula-difference evidence

## Authoritative state recovery

- Recovered run: `minimal6_states_live_20260713_v2`.
- Manifest SHA-256: `43305a3b4b3f565d600c2ccb783af36c62ed385e0ea6aa9c13f71b3a34e370c0`.
- State identity: `sketch/mirror x 0/-12/-19 V`; all six entries are passed and `outputs_complete=true`.
- Required raw fields: `ImpactIonization` (`cm^-3*s^-1`), `eVelocity`/`hVelocity` (`cm*s^-1`), and `eIonIntegral`/`hIonIntegral`/`MeanIonIntegral` (`1`).

The archive is read-only evidence. New exports record each state export's member SHA-256 table and reject altered members during validation.

A recovered-archive seal covering 355 members was generated and immediately
verified at
`build-release/reference_tcad/pn2d_sentaurus2018_minimal6/recovery_validation/minimal6_states_live_20260713_v2/recovery_validation.json`.
Its SHA-256 is
`c343539775337038848adc3a1f88b45880110d6960a3fe8ea3f3309032645f6d`,
and its recorded manifest hash is the value above.

## Scope

This is a diagnostic evidence package. It does not establish a physical BV curve or alter production solver defaults.
## Task 5 figure and documentation contract

The formula-difference CLI emits the fixed PNG/PDF pairs `gradient`, `current_alpha`, `source_waterfall`, `interaction`, and `topology_symmetry`, plus `figure_manifest.json`. The manifest retains units, the diagnostic disclaimer, and a reviewer/date/checklist entry. A visual inspection remains a Phase A gate; an automatically produced file is not evidence of a physical BV curve.

The reports maintain these source labels without substitution: native Sentaurus `ImpactIonization`, the separately labelled Sentaurus `alpha*|J|/q` reconstruction, and the separately labelled Vela `alpha*flux*partial_volume` reconstruction. Geometric zeros are marked as zeros rather than replacing them with log-floor values.

## Root-cause implementation map

| Classified factor or control | Sentaurus source | Independent diagnostic implementation | Production control implementation |
| --- | --- | --- | --- |
| `ni_eff/BGN`, alpha-law control | [models.par](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par), [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) | [physics.py](../../scripts/pn2d_minimal6_diagnostics/physics.py) | [ImpactIonizationModel.cpp](../../src/physics/ImpactIonizationModel.cpp) |
| Current semantics, gradient recovery, mobility, source mapping | [state deck](../../reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd) | [counterfactual.py](../../scripts/pn2d_minimal6_diagnostics/counterfactual.py), [support.py](../../scripts/pn2d_minimal6_diagnostics/support.py) | [DCSweep.cpp](../../src/simulation/DCSweep.cpp) |

Parameter agreement is a control, not a causal conclusion. The current formula report remains `insufficient_data` until raw Vela state/operator inputs permit declared counterfactual substitutions and closure.

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
