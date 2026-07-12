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

## 2026-07-11 GSS 2x2 Replay

The current-HEAD replay now covers two junction-doping strategies and two
avalanche-current discretizations:

- `legacy_density_gradient`
- `legacy_gss_midpoint`
- `reported_density_gradient`
- `reported_gss_midpoint`

The coarse Sentaurus TDR reports zero signed aggregate doping at the three
compensated junction nodes. The existing `dominant_signed_region` importer
therefore preserves those nodes by design. To reproduce the historical
`p -> p -> n` Vela baseline without changing the global importer, only the
isolated legacy variant copies replay unresolved `signed_aggregate_zero`
junction nodes as p-side nodes. The reported variants retain the raw
`p -> compensated -> n` import. The orchestrator rejects a run unless all four
meshes match, both variants within each doping strategy match, and the two
doping-strategy hashes differ.

Audited output:

`build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/pn2d_bv_compensated_gss_matrix_v2_20260711`

- Mesh SHA-256: `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e`.
- Legacy doping SHA-256: `af212731493ff2fb49be9224a9cecfe0f527561c906cc6b3ec01992f896233df`.
- Reported doping SHA-256: `714bb5c461d0acba49b1f9211318cc120a2e3891f92367776b636eda4b7fd155`.
- All four sweeps converged for 401 points through `-20 V`, with final
  handoff stage `newton`.
- The detailed matrix contains 72 rows and 36 matched
  `gss_midpoint / density_gradient` pairs.
- At `-12`, `-19`, and `-20 V`, GSS midpoint source ratios are approximately
  `0.995-0.996`, `0.995`, and `0.992-0.993`. The change follows the electron
  flux proxy; electron alpha and mobility ratios remain approximately one.
- The dominant residual classification remains
  `sg_discretization_ni_or_current_semantics`; midpoint reconstruction is a
  small correction, not a complete explanation of the Sentaurus/Vela gap.

## 2026-07-12 Triangle GSS GradQf Implementation

The opt-in source_mapping_mode triangle_gss_gradqf_truncated is implemented
with current_approximation cell_reconstructed, midpoint gss_logistic,
cell-gradient quasi-Fermi discretization, Genius-truncated geometry, and
symmetric endpoint partition.

For every Tri3 cell the implementation reconstructs two-dimensional
grad(phin) and grad(phip), evaluates alpha once per carrier and cell, evaluates
each local-edge current as mu * carrier_mid_gss * |Delta(phiF)| / h, multiplies
by the non-negative Genius-truncated partial area, and splits the result
equally between edge endpoints.

The midpoint orientation follows GSS 0.47:

    n_mid = n_i aux2((psi_j-psi_i)/(2 Vt))
          + n_j aux2((psi_i-psi_j)/(2 Vt))
    p_mid = p_i aux2((psi_i-psi_j)/(2 Vt))
          + p_j aux2((psi_j-psi_i)/(2 Vt))

This is not a byte-for-byte reconstruction of old GSS circumcentric geometry.
The alpha and current formulas follow GSS, while source geometry remains the
existing non-negative Genius-truncated policy. Non-Tri3 cells and surface
mobility are rejected rather than silently approximated.

Residual assembly, carrier diagnostics, VTK AvalancheGeneration, and the new
triangle_gss_sources CSV share the same per-cell evaluator. The coupled
Jacobian differentiates that evaluator against the three psi, three phin, and
three phip values of each cell. The dedicated audit improved from a pre-fix
relative difference of 9.767e-5 to 2.265e-8, below the 5e-5 gate.

### 2x3 Replay Matrix

Manifest v3 defines six stable variants:

- legacy_density_gradient
- legacy_gss_midpoint
- legacy_triangle_gss_gradqf
- reported_density_gradient
- reported_gss_midpoint
- reported_triangle_gss_gradqf

A complete matrix contains 108 detailed rows and 72 identified pairs. Triangle
records are aggregated across adjacent cells by global edge and merged with
the historical SG decomposition. A zero Genius-truncated edge area remains a
valid zero-source record.

Prepare-only evidence is under
build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/
pn2d_bv_compensated_gss_matrix_v3_prepare_20260712:

- six decks;
- manifest schema/version v3/3;
- 16 command records;
- only expected doping and source-discretization fields differ.

Full-run status evidence is under
build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/
pn2d_bv_compensated_gss_matrix_v3_full_status_20260712:

- density-gradient and historical GSS-midpoint variants complete 401 points;
- both triangle variants produce 390 points;
- last stable triangle bias is -19.4 V and failed trial bias is approximately
  -19.43650073 V;
- failure is line_search_non_decrease with positive finite carrier densities;
- manifest records both failed variants, outputs_complete=false, and a nonzero
  exit status;
- no 108-row report is generated without a converged -20 V triangle state.

At -19 V, legacy triangle total source is about 1.7233e21 versus 1.0944e15
for historical GSS midpoint, a factor of about 1.57e6. Field and alpha maxima
are comparable; the change is concentrated in GSS logistic midpoint
density/current support. Do not replace the default model or add a quasi-Fermi
clamp based on this result.

References: local PDF Chapter 9 equations 9.131-9.137; GSS 0.47 archive
SHA-256 f7359ea50ab19b8701dc241f990b8e462a2106b15eaf27613502bdf3472ba59d;
local Genius revision 543da845; local Charon revision 7cc387.

Next gate: same-state, same-edge audit of the GSS logistic midpoint against
GSS 0.47 and Sentaurus carrier/current vectors, especially edge 50 near
-19 V. Keep density_gradient as default and gss_midpoint as the lower-risk
experimental baseline until the 1.57e6 factor is explained.
