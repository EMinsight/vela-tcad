# TransportModels five-bias spatial oracle comparison

Status: **complete**. Node mapping: **exact**.

Profiles follow the silicon normal from the Si/SiO2 interface to 20 nm depth at the source end, channel midpoint, and drain end.

| Vg (V) | Qn p95 (mV) | n p95 (dex) | φn p95 (mV) | Enormal p95 (dex) | μn p95 (dex) | |∇φn| p95 (dex) |
|---:|---:|---:|---:|---:|---:|---:|
| -0.20 | 49.4977 | 0.886364 | 100.197 | 0.380389 | 0.112164 | 0.12586 |
| -0.04 | 38.3298 | 1.04443 | 96.9089 | 0.376038 | 0.0988905 | 0.128496 |
| 0.12 | 13.5579 | 0.50258 | 40.7496 | 0.366942 | 0.104588 | 0.132982 |
| 0.28 | 5.60432 | 0.235479 | 16.2662 | 0.344286 | 0.164108 | 0.139887 |
| 1.00 | 8.27027 | 0.177923 | 17.6821 | 0.109738 | 0.186564 | 0.177782 |

## Profile localization

| Vg (V) | Worst Qn profile / p95 (mV) | Worst n profile / p95 (dex) | Worst φn profile / p95 (mV) |
|---:|---:|---:|---:|
| -0.20 | drain_end / 52.921 | drain_end / 1.00482 | drain_end / 106.498 |
| -0.04 | drain_end / 39.2889 | drain_end / 1.05495 | drain_end / 99.9747 |
| 0.12 | drain_end / 15.2611 | drain_end / 0.544633 | drain_end / 45.6808 |
| 0.28 | drain_end / 6.69565 | drain_end / 0.251104 | drain_end / 17.5465 |
| 1.00 | drain_end / 8.88491 | drain_end / 0.334496 | drain_end / 22.7162 |

Qn criterion (`p95 <= 20 mV`): **fail**.
Electron-density criterion (`p95 <= 0.2 dex`): **fail**.

Profile figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\idvg_spatial_oracle_2026-08-21\comparison\idvg_spatial_profiles_vg1p00.png`
Error heatmap: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\idvg_spatial_oracle_2026-08-21\comparison\idvg_spatial_error_heatmap.png`
