# PN2D minimal6 formula-difference evidence

## Authoritative state recovery

- Recovered run: `minimal6_states_live_20260713_v2`.
- Manifest SHA-256: `b44ad95d5df6d57383ba3d5b292818568e358d67f0fc0424ee72f95b673e8aaa`.
- State identity: `sketch/mirror x 0/-12/-19 V`; all six entries are passed and `outputs_complete=true`.
- Required raw fields: `ImpactIonization` (`cm^-3*s^-1`), `eVelocity`/`hVelocity` (`cm*s^-1`), and `eIonIntegral`/`hIonIntegral`/`MeanIonIntegral` (`1`).

The archive is read-only evidence. New exports record each state export's member SHA-256 table and reject altered members during validation.

A recovered-archive seal covering 355 members was generated and immediately verified at uild-release/reference_tcad/pn2d_sentaurus2018_minimal6/recovery_validation/minimal6_states_live_20260713_v2/recovery_validation.json; its recorded manifest hash is the value above.

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
