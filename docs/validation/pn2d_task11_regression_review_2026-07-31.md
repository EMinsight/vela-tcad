# PN2D Task 11 regression and review

Date: 2026-07-31

Typed outcome:

```text
task11_regression_passed_keep_opt_in_independent_review_pending
```

## Decision

The Task 11 implementation and regression gates pass. The
`element_edge_sg_gss_laux` candidate remains opt-in because the historical
independent scientific and code reviews approved the operator only as a
diagnostic with unchanged production defaults. They do not cover the new M0
and M2 self-consistent golden-parity evidence or authorize a production
default change.

No physical parameter, numerical tolerance, empirical scale, voltage shift,
or production default was changed during Task 11.

## Scope and starting state

- Branch: `codex-pn2d-minimal6-operator-audit`
- Starting commit:
  `58d6c1725c2c2e06a5a1438d13af9e02a4544dba`
- Release runner SHA-256:
  `3d6acc8310246b4d46cd7d764ab171bab8d3302f843a46e70d634bff537e25e4`
- Generated simulation outputs remain under `build-release/` and are not
  committed.
- The user-owned untracked `tmp/` directory remains untouched.

## Regression results

| Gate | Required | Observed | Result |
|---|---:|---:|---|
| focused avalanche/source/Jacobian/Tri3/DC/continuity/PN2D CTest | zero failures | 130/130 passed | pass |
| full Release CTest | zero failures | 506/506 passed | pass |
| ASCII source check | pass | passed | pass |
| PN2D config template validation | pass | passed | pass |
| Sentaurus import/schema validation | pass | passed | pass |
| BV process observability/manifest validation | pass | passed | pass |
| forward IV convergence | 201/201 | 201/201 | pass |
| maximum forward anchor change | `< 0.5%` | `0.00178944%` | pass |
| IV mobility doping basis | `cell_reconstructed_total_impurity` | exact match | pass |
| BV mobility doping basis | `net_doping` | exact match | pass |
| avalanche-off nonzero-bias RMSE | `<= 0.001 dex` | `6.9965091e-5 dex` | pass |
| avalanche-off maximum error | `<= 0.002 dex` | `1.2546889e-4 dex` | pass |
| generated output or `tmp/` staged | none | none | pass |
| new independent scientific and code review | both accept | not yet obtained | pending |

Commands:

```powershell
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure `
  -R "(avalanche|impact|GSS Laux|continuity|DCSweep|pn2d|ascii_sources)"
ctest --test-dir build-release --output-on-failure
ctest --test-dir build-release --output-on-failure `
  -R "^(ascii_sources|pn2d_config_templates|sentaurus_import_tools|pn2d_bv_process_observability_regression)$"
```

The first targeted validator invocation reached the local command timeout
after `ascii_sources` and `pn2d_config_templates` had passed. The remaining
`sentaurus_import_tools` and
`pn2d_bv_process_observability_regression` tests were rerun explicitly and
both passed. The complete 506-test run had already passed all four tests.

## Forward-IV confirmation

The 201-point `0..20 V` experiment was regenerated from the sealed
configuration:

```powershell
python scripts\run_pn2d_forward_mobility_doping_basis_experiment.py `
  --runner build-release\vela_example_runner.exe `
  --baseline-config build-release\pn2d-forward-mobility-doping-basis-20260727\cell_reconstructed_total_impurity\simulation_cell_reconstructed_total_impurity.json `
  --sentaurus-fields build-release\pn2d-forward-field-audit-20260727\sentaurus_fields `
  --out-dir build-release\pn2d-task11-forward-mobility-basis-20260731
```

All three observation candidates converged at 201/201 points. The contractual
candidate is `cell_reconstructed_total_impurity`.

| Bias (V) | Sealed current (A/um) | Task 11 current (A/um) | Change |
|---:|---:|---:|---:|
| 1 | `1.1476888366e-4` | `1.1476682994e-4` | `-0.00178944%` |
| 2 | `9.4193508862e-4` | `9.4192975146e-4` | `-0.000566617%` |
| 5 | `3.2809770786e-3` | `3.2809717898e-3` | `-0.000161197%` |
| 10 | `7.0866107001e-3` | `7.0866055698e-3` | `-0.0000723954%` |
| 15 | `1.0791901810e-2` | `1.0791896815e-2` | `-0.0000462909%` |
| 20 | `1.4403941547e-2` | `1.4403936676e-2` | `-0.0000338130%` |

The relative change at `0.1 V` is not used as an acceptance metric because
both currents are near the numerical zero-current floor. The predeclared
nonzero-current anchors above all pass.

The regenerated contractual candidate also retains:

- exact-anchor median absolute error versus Sentaurus: `0.271431%`;
- current error at `20 V` versus Sentaurus: `-0.262471%`;
- electron/hole junction-mobility RMSE:
  `8.69528/1.87221 cm^2/(V s)`.

Hashes:

```text
configuration:
fd0272b9f7ff72547a60d074ff5b92cc2b9fb9d80a9646c66ab758c891a14677
mesh:
c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e
doping:
714bb5c461d0acba49b1f9211318cc120a2e3891f92367776b636eda4b7fd155
materials:
212b896ed86ace76d0d02e86b90f6ced6a2851417b2ec321d0dd702cfbfa6524
Task 11 curve:
475a7405879b9c4539b7a9f6ada3c66d2bcec3d05a680ad7923a96bce0362cee
candidate summary:
e2c525414daf3f019dcffa54e6d7561480e4a3ce84e0d5b96af877563322c394
```

## BV/IV configuration contract

Fresh default configurations and manifests were rendered with
`scripts/generate_pn2d_config.py`.

| Template | Mobility model | Doping basis | Impact model | Sweep |
|---|---|---|---|---|
| `pn2d_iv` | `masetti` | `cell_reconstructed_total_impurity` | `none` | `iv` |
| `pn2d_bv` | `masetti_field` | `net_doping` | `van_overstraeten` | `bv_reverse` |

The production BV template still selects
`current_approximation=cell_reconstructed`. The SG/Laux candidate is enabled
only by the exact-lattice experiment runner, which explicitly writes
`current_approximation=element_edge_sg_gss_laux`.

Generated contract hashes:

```text
pn2d_iv.json:
26abf7715699a73c11799e73cfd9a376e7790d7f203861aa0405b73bbcc8c640
pn2d_iv.manifest.json:
432502095ac353585688ad5bfc3aa4639a1595ab477002d586161d57f5bc2899
pn2d_bv.json:
ecedb41106879209cb241369f5edbcfef9e867621203bd58e2502a4ff28bbc22
pn2d_bv.manifest.json:
b9c15381bcbd7065da0fc33030920c7ab37bc8daa4630906c4810cc5fea314cf
```

## Avalanche-off M0 confirmation

The current Release runner regenerated the sealed 21-point M0
avalanche-off/SRH audit against the original paired Sentaurus curve and TDR
sequence:

```powershell
python scripts\run_pn2d_bv_off_srh_spatial_audit.py `
  --out-dir build-release\pn2d-task11-bv-off-srh-paired-20260731 `
  --sentaurus-curve build-release\pn2d-bv-off-srh-mesh-matrix-coarse-baseline-20260728\paired_audits\M0\sentaurus_curve.csv `
  --sentaurus-source build-release\pn2d-bv-off-srh-mesh-matrix-coarse-baseline-20260728\sentaurus_runs\pn2d_srh_coarse_mesh_M0_20260728\source `
  --sentaurus-intervals 400
```

Results:

- Vela/Sentaurus coverage: `21/21` and `21/21`;
- 20-point nonzero-bias log-current RMSE: `6.9965091e-5 dex`;
- maximum log-current error: `1.2546889e-4 dex`;
- maximum electron/hole source-contact closure:
  `7.4819873e-7/3.5611877e-9`;
- maximum terminal-pair closure: `3.7621206e-23 A/um`;
- current Vela curve SHA-256 equals the sealed 2026-07-28 curve SHA-256:
  `18a12f632813e0e826f75ce91796398de479d30e9b35a9409be8b9d2c7e7335f`.

Run and report manifest hashes:

```text
run_manifest.json:
f085b00aae3e41f11f695379bb4b4fb2cce3525dca31d2c0fdf854716eebc1ac
report/report_manifest.json:
291fdcd488c48b47fb09d0aa0ebd4c8ede4c582b79de3a77dd8cfe0b44621768
```

An initial Task 11 invocation used the script's current default
`sentaurus_on_off_refresh_integer.csv`. Its hash and current levels differ
from the sealed Task 9 M0 reference, so that comparison was rejected before
acceptance scoring. The valid rerun above binds the exact sealed curve hash
and the matching 400-interval TDR sequence.

## Scientific evidence review

The current evidence supports the following limited claims:

1. The production avalanche coefficients, QFP-gradient drive, mobility
   parameters, and physical inputs were not fitted.
2. SG/Laux changes the element-local current/source support semantics; it does
   not change the Van Overstraeten coefficient law.
3. The stable same-grid M0 and M2 comparisons support Vela-versus-Sentaurus
   golden parity. M0 has `|delta V_break|=0.032 V`; M2 has
   `|delta V_break|=0.014 V`, median current error `0.05663 dex`, and maximum
   current error `0.07953 dex`.
4. The shared M1 nonmonotonic branch is retained as a topology/continuation
   observation. It is not evidence that Vela alone differs from Sentaurus.
5. M1-to-M2 source and knee changes still mean that a mesh-converged physical
   breakdown knee has not been demonstrated. Same-grid golden parity and
   cross-grid physical convergence are separate claims.
6. The low-current `-3..-7 V` nonmonotonicity remains classified as a
   residual-floor/state-precision effect rather than an avalanche-source or
   terminal-current-extraction defect.

This review finds no regression or contradiction requiring rejection of the
opt-in candidate. It does not replace an independent review of the new
self-consistent evidence.

## Code evidence review

- Task 11 made no solver-source change.
- The SG/Laux branch remains explicitly selected and guarded by its canonical
  element-box configuration checks.
- The production C++ default remains `mobility_density_gradient`; the
  production PN2D BV template remains `cell_reconstructed`.
- BV and IV doping bases remain separated and are protected by generated
  configuration tests.
- The Release test suite, focused physics tests, ASCII check, schema/import
  tests, process-observability tests, and forward/off numerical baselines all
  pass.
- No generated simulation artifact is included in the source diff.

This review finds no code regression and no accidental default change. A
fresh independent code reviewer must still approve any later proposal to
change the production default.

## Independent-review boundary

The 2026-07-25 historical independent scientific and code reviews approved
the SG/Laux implementation with the explicit conclusion
`production default unchanged`. They did not evaluate the later
self-consistent M0/M2 parity, the corrected balanced-junction contract, or a
new default proposal.

Therefore the Task 11 regression work is complete, but the final acceptance
criterion requiring two new independent approvals is intentionally left
pending. The allowed decision is:

```text
keep element_edge_sg_gss_laux opt-in;
do not change the production default;
preserve the M1 mesh/continuation observation;
request separate scientific and code review before any default-change patch.
```
