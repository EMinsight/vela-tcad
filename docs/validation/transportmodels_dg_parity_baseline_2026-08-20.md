# TransportModels DG parity-improvement frozen baseline

Date: 2026-08-20

Status: **frozen and reproducible**. This baseline is the phase-0 control for
improving Vela/Sentaurus DG parity. It does not claim that the improvement
targets already pass.

## Frozen curve metrics

| Metric | Current value | Improvement target |
|---|---:|---:|
| DD Id-Vd maximum relative error | 0.866356% | <= 2.00% |
| DD Id-Vg on-state maximum relative error | 18.495522% | <= 20.00% |
| DG Id-Vd maximum relative error | 10.721943% | <= 5.00% |
| DG Id-Vd endpoint relative error | 9.842376% | <= 3.00% |
| DG Id-Vg on-state maximum relative error | 20.785387% | <= 10.00% |
| DG Id-Vg transition maximum log error | 0.339151 dex | <= 0.150 dex |
| DG Qn surface 99th-percentile absolute error | 186.867294 mV | phase-1 target uses p95 <= 20.0 mV |
| DG electron-density surface 99th-percentile log error | 2.418067 dex | phase-1 target uses p95 <= 0.20 dex |

## Reproduction

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\freeze_transportmodels_dg_parity_baseline.py --check
```

The JSON companion records absolute paths and SHA-256 values for every frozen
candidate, final state, workflow manifest, regional analysis, and spatial
summary. Curve files must contain 21 unique finite points aligned to the
Sentaurus bias lattice within 1e-12 V.

## Policy

- DD remains a non-regression control while DG is modified.
- Deep-off Id-Vg is excluded from ordinary relative-error gates.
- No later phase may relax nonlinear convergence tolerances to pass an
  amplitude target.
- A changed artifact hash requires an explicitly regenerated phase-0 baseline,
  not silent acceptance by `--check`.
