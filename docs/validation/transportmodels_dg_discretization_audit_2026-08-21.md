# TransportModels DG discretization audit

The audit evaluates all supported global DG operators on the same fixed state,
corrected material contract, and neutral interface. Raw residual magnitudes are
reported for hotspot localization, but cannot by themselves rank formulations
with different primary variables or row scaling. A self-consistent current test is
therefore required before selecting the production operator.

| Discretization | Diagnostic | Global raw L1 | Max raw residual | Max node |
|---|---|---:|---:|---:|
| P1 direct | diagnostic_complete | 483220 | 21372.4 | 2 |
| P1 lambda | diagnostic_complete | 140514 | 7376.3 | 69 |
| CVFEM full | diagnostic_complete | 483220 | 21372.4 | 2 |
| Sentaurus box | diagnostic_complete | 405275 | 3391.67 | 1640 |
| GSS potential-like | diagnostic_complete | 564137 | 24744.5 | 2 |
| GSS density | diagnostic_complete | 130907 | 7331.7 | 123 |
| Conservative sqrt | diagnostic_complete | 2.1668e+21 | 6.02462e+20 | 1703 |

## Decision

- `p1_direct` remains the conservative control because the phase-2/3 residual
  decomposition and its units are already audited.
- `sentaurus_box` with a neutral interface advances as the primary contender:
  it preserves the potential-like variable and substantially reduces both audited
  fixed-state metrics relative to corrected `p1_direct`.
- Fitted density/square-root formulations remain diagnostic candidates until their
  self-consistent convergence and terminal current are demonstrated.

Figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_discretization_fixed_state_audit_2026-08-21\discretization_residuals.png`
CSV: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_discretization_fixed_state_audit_2026-08-21\discretization_residuals.csv`
