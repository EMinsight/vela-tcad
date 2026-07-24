# PN2D Minimal6 high-field box-current replay

## Result

The typed `sentaurus_lowfield_element_electric_field` branch passed the
edge-current, sign, and terminal-current gates for all 40 exact states.  It
retained a bounded fixed-state KCL failure, so the result is classified as
`bounded_gate_failure`, not as a self-consistent solution.

| Quantity | Electron | Hole |
|---|---:|---:|
| valid active edges | 200 | 200 |
| median absolute current error (dex) | 0.000853883 | 0.000446990 |
| P95 absolute current error (dex) | 0.001343436 | 0.000570236 |
| maximum absolute current error (dex) | 0.001348167 | 0.000576072 |
| sign agreement | 100% | 100% |

The maximum reconstructed anode-current relative error is
`1.428545e-3` (0.143%).  The maximum internal-node current divergence
relative to terminal current is `8.731080e-4`; the strict `1e-8` KCL gate
therefore fails.  This is the expected diagnostic consequence of changing
mobility while holding the imported Sentaurus state fixed.

## Determinism and independent verification

The two roots

- `build-release/pn2d-minimal6-highfield-box-current-20260724-a`
- `build-release/pn2d-minimal6-highfield-box-current-20260724-b`

are byte-identical for every generated CSV, report, and independent
verification JSON.  Both independent verifications report
`passed_expected_bounded_gate_failure` with zero failures.

No production formula was changed in this task.
