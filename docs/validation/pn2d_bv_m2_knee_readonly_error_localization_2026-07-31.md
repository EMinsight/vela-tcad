# PN2D M2 BV knee-region read-only error localization

## Outcome

The growing `avalanche_on` discrepancy from -18 V through -20 V is localized to the self-consistent carrier-density and SG/Laux current-amplitude loop. It is not localized to the ionization coefficient, QFP-gradient driving force, dominant-edge mobility, carrier partition, source geometry measure, or hotspot movement.

This is an observation-only result. No solver setting, physical model, PN2D template default, or acceptance threshold was changed.

## Quantitative evidence

| Quantity, -18 V to -20 V | Change |
|---|---:|
| Terminal-current absolute error | +0.05521 dex |
| Integrated-source absolute error | +0.05488 dex |
| Electron active-cell current deficit | +0.05319 dex |
| Hole active-cell current deficit | +0.05314 dex |
| Electron source-weighted density deficit | +0.05135 dex |
| Hole source-weighted density deficit | +0.05608 dex |

The bias-wise correlation between terminal-current error and integrated-source error is `0.99999235`. Correlations with the electron and hole source-weighted density errors are `0.99790880` and `0.99998746`.

The total source deficit grows from `0.02866 dex` at -18 V to `0.08354 dex` at -20 V. The terminal-current error grows from `0.02361 dex` to `0.07883 dex` over the same points. Electron and hole source fractions remain paired within `0.00263`, so the discrepancy is not carrier-specific.

## Controls that do not explain the growing error

- QFP-gradient drive: maximum active-cell mismatch `0.00177 dex`.
- Ionization coefficient: source-weighted integral counterfactual within `0.00854 dex`.
- Dominant-edge mobility: maximum active-cell mismatch `0.04556 dex`; it does not reproduce the common electron/hole error growth.
- Source measure: Sentaurus and Vela element-vertex measures agree within `9.47e-6` relative on every nonzero compared support.
- Hotspot: the maximum-source cell and element vertex are identical at all 11 knee biases.
- Spatial mapping: normalized overlap is at least `0.836` at element-vertex support and `0.902` after cell aggregation. Overlap improves as current error worsens, excluding a growing geometric remap or hotspot relocation as the primary driver.
- Determinism: `iv.csv`, `process_probe.csv`, `newton_attempts.csv`, and `newton_history.csv` are byte-identical between independent Vela runs A and B.

## Edge-current interpretation

Edges were paired by unordered endpoint connectivity because Sentaurus and Vela use different local-edge numbering for the same triangle.

The median Vela-versus-native-Sentaurus edge-current projection error grows from about `0.18 dex` at -18 V to `0.24 dex` at -20 V, tracking the knee discrepancy. The Sentaurus operator-replay field has a stable `1e6` scale signature: its raw comparison is about `6 dex` off, whereas the inferred square-micrometre to square-centimetre conversion reduces the median mismatch to `0.10-0.12 dex`. Both raw and corrected values are retained; the correction is not applied to production physics.

## Interpretation and limitation

The evidence locates the growing error in the self-consistent carrier/current/source feedback path. Drive, alpha, geometry, and carrier partition are stable controls. However, paired self-consistent snapshots cannot establish which side of the feedback loop moves first: the carrier-density deficit or the native-current versus SG/Laux support difference.

The next discriminating experiment should remain read-only: import the exact M2 Sentaurus states at -18, -19.5, -19.7, and -20 V, run Vela SG/Laux with `coupling_mode=postprocess_only`, and compare its integrated source with both Sentaurus and self-consistent Vela. No acceptance threshold should be changed.

## Reproduction

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\diagnose_pn2d_bv_m2_knee_error.py `
  --vela-manifest build-release\pn2d-bv-template-default-prospective-v2-default-20260731\M2\run-a\manifest.json `
  --vela-probe build-release\pn2d-bv-template-default-prospective-v2-default-20260731\M2\run-a\avalanche_on\process_probe.csv `
  --vela-run-a-root build-release\pn2d-bv-template-default-prospective-v2-default-20260731\M2\run-a `
  --vela-run-b-root build-release\pn2d-bv-template-default-prospective-v2-default-20260731\M2\run-b `
  --sentaurus-manifest build-release\pn2d-task10-balanced-m2-sentaurus-process-v2-20260731\manifest.json `
  --parity build-release\pn2d-bv-template-default-prospective-v2-default-20260731\M2\contract-domain-parity.json `
  --output-root build-release\pn2d-bv-m2-knee-readonly-diagnostic-20260731

D:\msys64\ucrt64\bin\python.exe scripts\verify_pn2d_bv_m2_knee_error.py `
  --root build-release\pn2d-bv-m2-knee-readonly-diagnostic-20260731
```

Machine-readable outputs are written below `build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/` and remain generated, untracked artifacts.
