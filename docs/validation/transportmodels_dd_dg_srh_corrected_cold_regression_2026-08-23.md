# TransportModels corrected cold-start DD/DG regression

Execution status: **complete**; acceptance status: **fail**; completed `84/84` comparison points.

Sweep semantics match the Sentaurus decks: Id-Vg initializes at `Vg=-1 V`, ramps the drain to `1.1 V`, then sweeps to `2.2 V`; Id-Vd initializes separately at `Vg=1 V`, `Vd=0 V`, then sweeps to `2 V`.

| Branch | Curve | Points | Primary metric | Secondary metric | Tertiary metric | Numerical qualification |
|---|---|---:|---:|---:|---:|---|
| DD | idvg | 21/21 | off 0.0677308 dex | transition 0.0318431 dex | on 2.9281% | unresolved at [-1.0, -0.84, -0.68] |
| DD | idvd | 21/21 | max 1.3613% | median 1.3449% | endpoint 1.3429% | resolved |
| DG | idvg | 21/21 | off 0.316766 dex | transition 0.43823 dex | on 9.1535% | unresolved at [-1.0, -0.84, -0.68] |
| DG | idvd | 21/21 | max 8.6692% | median 7.9974% | endpoint 7.7250% | resolved |

Figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_dg_srh_corrected_cold_regression_2026-08-23\dd_dg_idvg_idvd_comparison.png`
