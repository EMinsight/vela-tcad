# TransportModels DG validated quantum-contract regression

Execution: **complete**; main-curve acceptance: **pass**; completed `42/42` points.

The regression uses the corrected material/Fermi-BGN/SRH contract and the independently validated `include_insulators + sentaurus_box` DG contract. Deep-off Id-Vg is reported separately and does not veto the normal-region result.

| Curve/region | Metric | Result | Limit | Status |
|---|---|---:|---:|---|
| Id-Vg transition | max log error | 0.0328737 dex | 0.15 dex | pass |
| Id-Vg on | max relative error | 2.710% | 10.0% | pass |
| Id-Vd | max relative error | 2.432% | 5.0% | pass |
| Id-Vd 2 V | endpoint relative error | 0.976% | 3.0% | pass |

## Improvement from the corrected cold baseline

| Metric | Prior DG baseline | Validated quantum contract |
|---|---:|---:|
| Id-Vg transition max log error | 0.43823 dex | 0.0328737 dex |
| Id-Vg on max relative error | 9.154% | 2.710% |
| Id-Vd max relative error | 8.669% | 2.432% |
| Id-Vd 2 V endpoint error | 7.725% | 0.976% |

## Separate deep-off branch

| Vg (V) | Log error (dex) | Id/KCL residual | Classification |
|---:|---:|---:|---|
| -1.00 | 0.00426475 | 378.273 | pass |
| -0.84 | 0.00476344 | 377.241 | pass |
| -0.68 | 0.0115386 | 537.642 | pass |

## Id-Vd continuation history

The extra points below are solver-path bridges only; they are excluded from the strict 21-point comparison lattice.

| Restart Vd (V) | Bridge Vd (V) | Bridge converged | Attempt final status |
|---:|---:|---|---|

Figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_srh_density_coupling_sentaurus_default_2026-08-23\full_21_point_curves\dg_quantum_contract_idvg_idvd_comparison.png`

Run manifest: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_srh_density_coupling_sentaurus_default_2026-08-23\full_21_point_curves\runs\dg\workflow_manifest.json`
