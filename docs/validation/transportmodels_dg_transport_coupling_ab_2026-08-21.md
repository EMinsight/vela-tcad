# TransportModels DG transport-coupling A/B

Status: **partial**.

| Vg (V) | Vela exponential/direct (dex) | Sentaurus default/DirectQC (dex) | residual effect error (dex) |
|---:|---:|---:|---:|
| -0.20 | 0.995653 | 0.005193 | 0.99046 |
| -0.04 | 0.927603 | 0.00453061 | 0.923072 |
| 0.12 | failed | 0.0035936 | failed |
| 0.28 | 0.260495 | 0.00118659 | 0.259308 |
| 1.00 | 0.0066905 | -0.0040827 | 0.0107732 |

The implementation intentionally leaves `direct_band_edge` as the compatibility default.
