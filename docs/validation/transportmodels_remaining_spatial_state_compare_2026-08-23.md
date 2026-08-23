# TransportModels remaining DG spatial-state comparison

Corrected self-consistent Vela states are compared node-for-node with the new Sentaurus 2022 TDR snapshots.

| Group | Bias | Qn p95 all/drain (mV) | n p95 all/drain (dex) | phin p95 all/drain (mV) |
|---|---:|---:|---:|---:|
| dg_idvd | 0.20 V | 142.4 / 159.1 | 0.4603 / 1.007 | 1.067 / 1.439 |
| dg_idvd | 0.50 V | 140.1 / 147 | 0.4715 / 1.194 | 6.292 / 2.247 |
| dg_idvd | 1.00 V | 137.9 / 137.4 | 0.5917 / 1.294 | 7.892 / 24.64 |
| dg_idvd | 2.00 V | 138.4 / 118.3 | 0.8516 / 1.418 | 16.67 / 48.76 |
| dg_idvg_deep_off | -1.00 V | 141 / 445.9 | 0.9059 / 5.053 | 8.711 / 82.7 |
| dg_idvg_deep_off | -0.84 V | 136.7 / 350.5 | 0.8095 / 4.689 | 8.051 / 46.53 |
| dg_idvg_deep_off | -0.68 V | 134.5 / 321.8 | 0.7837 / 4.245 | 8.186 / 47.67 |

Raw report: `D:\code-repo\vela-tcad\docs\validation\transportmodels_remaining_spatial_state_compare_2026-08-23.json`
