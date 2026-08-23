# TransportModels DG surface-mobility audit

Frozen Sentaurus quantum potential at `Vg=1 V`, `Vd=2 V`; all cases use the
corrected material contract. This isolates the classical mobility/transport path.

| Variant | Converged | Id error | Channel median Enormal (V/cm) | Mobility (cm2/Vs) | Drive (V/cm) |
|---|---:|---:|---:|---:|---:|
| Lombardi implicit interfaces | True | 1.3156% | 250642 | 18.9381 | 8893.44 |
| Lombardi explicit channel | True | 2.0399% | 250636 | 18.9379 | 8939.14 |
| No Enormal | True | 17.9055% | 250419 | 47.0582 | 10204.3 |
| No high-field saturation | True | 170.7999% | 248244 | 19.0712 | 22950.4 |

## Decision

- Best Frozen-Q terminal-current agreement: **Lombardi implicit interfaces**.
- The explicit channel selector is preferred only if it improves current and spatial
  mobility agreement; otherwise the existing model remains frozen to avoid fitting
  one endpoint at the expense of the curve.
- Vela Enormal is reconstructed from the adjacent substrate triangle and the exact
  R.Substrate/R.Gateox interface normal; Sentaurus values are native `eEnormal`.

Figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_surface_mobility_frozen_q_2026-08-21\surface_mobility_comparison.png`
