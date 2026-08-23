# TransportModels Sentaurus fixed-state formula replay

Sentaurus 2022 fields are inserted into Vela production density, mobility, SG flux, and SRH operators without a nonlinear update.

| Group | Bias | n p95 (dex) | mobility p95 (all/active) | SG J p95 (all/active) | SRH shape TV |
|---|---:|---:|---:|---:|---:|
| dg_idvg_transition | -0.20 V | 0.008396 | 0.4052 / 0.3293 | 0.767 / 0.5494 | 0.09081 |
| dg_idvg_transition | -0.04 V | 0.008357 | 0.3864 / 0.2756 | 0.8396 / 0.4662 | 0.08871 |
| dg_idvg_transition | 0.12 V | 0.008352 | 0.3688 / 0.2328 | 0.872 / 0.4379 | 0.08776 |
| dg_idvg_transition | 0.28 V | 0.008537 | 0.3644 / 0.235 | 0.8867 / 0.426 | 0.08699 |
| dg_idvg_transition | 1.00 V | 0.008129 | 0.3446 / 0.2239 | 0.8945 / 0.3584 | 0.08676 |
| dg_idvd | 0.20 V | 0.007971 | 0.303 / 0.2183 | 0.7885 / 0.2944 | 0.09011 |
| dg_idvd | 0.50 V | 0.008032 | 0.3331 / 0.1769 | 0.8547 / 0.3055 | 0.08856 |
| dg_idvd | 1.00 V | 0.008117 | 0.3442 / 0.2134 | 0.8858 / 0.3451 | 0.08568 |
| dg_idvd | 2.00 V | 0.008267 | 0.3516 / 0.2511 | 0.904 / 0.4528 | 0.08214 |
| dg_idvg_deep_off | -1.00 V | 0.008697 | 0.4586 / 0.3886 | 1.97 / 1.617 | 0.1354 |
| dg_idvg_deep_off | -0.84 V | 0.008665 | 0.4572 / 0.4087 | 1.875 / 1.544 | 0.1235 |
| dg_idvg_deep_off | -0.68 V | 0.008588 | 0.4443 / 0.392 | 1.43 / 1.42 | 0.1129 |

The SG comparison projects Sentaurus nodal current density onto Vela primal edges and integrates over the Vela dual couple. It is a localization diagnostic rather than an assertion that the two discretizations are identical.

Raw artifact directory: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_sentaurus_formula_replay_20260823`
