# PN2D BV Compensated Junction Source Proxy Debug Plan

Date: 2026-07-09

## Summary

The next BV debug step is evidence collection, not solver modification. The immediate goal is to explain why the compensated-junction probe removes the left/right `phin` drop asymmetry but still leaves a right-heavy avalanche source near the junction.

Known state:

- Baseline Vela coarse junction representation: `p -> p -> n`.
- Sentaurus coarse junction representation: `p -> compensated -> n`.
- Compensated probe representation: `p -> compensated -> n` on the Vela x=1.0 um junction column.
- Compensated probe changes median right/left `phin` drop ratios from `120x / 24x / 17x` to `0.8906 / 0.9291 / 0.9323` at `-12 / -19 / -20 V`.
- Remaining source right/left ratios are `17.95 / 5.80 / 4.89`, so there is still a source-proxy problem after QF-drop balancing.

## Diagnostic Added

Script:

```powershell
scripts\diagnose_pn2d_bv_compensated_source_proxy.py
```

Run command:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe scripts\diagnose_pn2d_bv_compensated_source_proxy.py `
  --baseline-report-root build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\bv_density_gradient_aligned_20260709 `
  --probe-root build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\bv_density_gradient_aligned_20260709\compensated_junction_probe_20260709 `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\bv_density_gradient_aligned_20260709\sentaurus_multibias_coarse_0p05 `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\bv_density_gradient_aligned_20260709\compensated_junction_proxy_compare_20260709
```

Outputs:

- `compensated_source_proxy_compare.csv`
- `compensated_source_proxy_compare_summary.json`
- `compensated_source_proxy_compare_report_20260709.md`

Expected coverage:

- `3 bias points x 3 y-cuts x 2 sides x 2 variants = 36 rows`
- Bias points: `-12 V`, `-19 V`, `-20 V`
- Y cuts: `0`, `0.25`, `0.5 um`
- Variants: baseline density-gradient, compensated-junction probe

## Current Finding

The compensated probe balances `phin` drops, but it does not remove the source right bias.

Median right/left ratios after compensation:

| bias | total source | electron source | hole source | electron alpha | electron flux proxy | electron alpha x flux | electron mobility |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -12 | 17.95 | 27.22 | 0.0546 | 0.449 | 60.64 | 27.22 | 1.117 |
| -19 | 5.80 | 21.71 | 0.212 | 0.727 | 29.88 | 21.71 | 1.074 |
| -20 | 4.89 | 17.85 | 0.246 | 0.749 | 23.84 | 17.85 | 1.071 |

Interpretation:

- The residual right-heavy source is carried by the electron source channel.
- The hole source channel is left-heavy at all three inspected biases.
- Electron alpha is below 1 right/left, and electron mobility is near unity.
- The remaining right-heavy multiplier follows electron SG flux proxy / raw flux proxy.
- Therefore the next debug target is density-gradient SG current/source construction and carrier-density/flux proxy selection, not a QF hard limiter.

## Policy On Doping Classification

`donors - acceptors` may be used to classify nodes and edges as `p`, `n`, `compensated`, `p-p`, `p-compensated`, `compensated-n`, or `p-n` for diagnostics, artifact alignment, and report labeling.

Do not directly clamp, hard-zero, or truncate `phin/phip` based on this classification. A solver-level limiter should only be considered after independent failure evidence and a targeted regression test show that it is necessary.

Default compensated threshold:

```text
abs(donors - acceptors) <= 1e-6 * max(abs(donors), abs(acceptors), 1.0)
```

## Next Debug Steps

1. Replay the density-gradient SG source construction at `-12 V`, `-19 V`, and `-20 V` for the same left/right junction edges.
2. Split electron flux proxy into endpoint densities, midpoint/reconstructed density, Bernoulli factor, raw signed flux, final flux proxy, and source weighting.
3. Compare `p-compensated` vs `compensated-n` edges for carrier-density support and SG exponential support.
4. Check whether Vela's density-gradient reconstruction creates an electron-density support mismatch relative to the Sentaurus coarse compensated junction.
5. Only after the SG flux/source evidence closes should solver-level branch limiters or QF clamping be discussed.

## Validation

Required checks for this diagnostic/documentation step:

- The script runs successfully and writes 36 CSV rows.
- The markdown report exists and states the dominant residual source factor.
- The handoff document includes the no-clamp/no-zero rule for `phin/phip`.
- `git diff --check` passes.
## Execution Status

Recorded on 2026-07-10 before local branch commit:

- Diagnostic script completed successfully and wrote `36` CSV rows.
- `compensated_source_proxy_compare_report_20260709.md` exists and records the electron SG flux proxy / raw flux proxy as the remaining right-heavy source factor after QF-drop balancing.
- The handoff document records the no-clamp/no-zero rule for `phin/phip`.
- `python -m py_compile scripts/diagnose_pn2d_bv_compensated_source_proxy.py` passed.
- `git diff --check` passed.
