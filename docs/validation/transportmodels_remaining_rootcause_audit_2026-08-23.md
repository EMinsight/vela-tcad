# TransportModels remaining-error root-cause audit

## Scope

This audit covers the three failures left after the corrected 84-point DD/DG
regression: DG Id-Vg transition behavior, systematically high DG Id-Vd, and
the first three deep-off Id-Vg points. Sentaurus Device T-2022.03-SP2 was used
to export seven additional 3315-node TDR states:

- DG Id-Vd at Vd = 0.2, 0.5, 1.0, and 2.0 V with Vg = 1.0 V.
- DG Id-Vg at Vg = -1.0, -0.84, and -0.68 V with Vd = 1.1 V.

Together with the existing five DG transition states, twelve Sentaurus states
were replayed through Vela production density, mobility, SG-flux, and SRH
operators without a nonlinear update.

## Fixed-Sentaurus-state formula replay

| Evidence | Result |
|---|---:|
| Electron-density reconstruction p95 | 0.00797-0.00870 dex |
| Hole-density reconstruction p95 | about 0.0393 dex |
| Electron-mobility median error | about 0.016 dex |
| Electron-mobility p95 error | 0.30-0.46 dex |
| Normal/transition current-carrying SG median error | about 0.03 dex |
| Deep-off current-carrying SG median error | 0.21-0.27 dex |
| SRH normalized shape total-variation distance | 0.082-0.135 |
| SRH top-50 hotspot overlap | 84%-94% |

Replacing Vela's reconstructed n/p by the exact Sentaurus n/p at the three
deep-off points changes the Vela SRH integral by only about 0.1% and changes
the carrier residual norm by less than about 1% in the resolved cases. The
density/Fermi/BGN reconstruction is therefore not the leading remaining error.

The mobility error is localized. Representative high-error current-carrying
edges lie at the silicon surface and source/drain extensions. At Vg = 1 V,
Vd = 2 V, several drain-surface edges have Sentaurus electron mobility of
roughly 11-32 cm2/(V s), while Vela gives about 0.5-2.7 cm2/(V s). Because the
baseline Vela Id is nevertheless high, excessive mobility is not the cause of
the high Id-Vd curve.

The SG edge diagnostic initially exposed a fixed two-decade offset. This was
a diagnostic-only unit bug: the native line flux was labeled as particles per
metre without applying the same line-integral factor used by terminal current.
The output conversion and a focused unit test were added; solver currents and
the 84-point curves were not changed.

## Self-consistent state comparison

With the current 84-point DG configuration, the Id-Vd Sentaurus/Vela state
difference grows with drain bias:

| Vd | Qn p95 | n p95 | drain-profile n p95 | drain-profile phin p95 |
|---:|---:|---:|---:|---:|
| 0.2 V | 142.4 mV | 0.460 dex | 1.007 dex | 1.44 mV |
| 0.5 V | 140.1 mV | 0.472 dex | 1.194 dex | 2.25 mV |
| 1.0 V | 137.9 mV | 0.592 dex | 1.294 dex | 24.6 mV |
| 2.0 V | 138.4 mV | 0.852 dex | 1.418 dex | 48.8 mV |

At the three deep-off points, the drain-profile electron-density p95 error is
4.24-5.05 dex and Qn p95 is 322-446 mV. These large local state errors coexist
with terminal currents near the cancellation/noise floor.

## Quantum-contract A/B

The earlier spatial baseline used `include_insulators=true`,
`global_discretization=sentaurus_box`, and the calibrated coefficient mass
ratio. The current 84-point regression had reverted to
`include_insulators=false` and omitted the Sentaurus-box discretization.

Restoring only that DG contract while keeping the corrected material, Fermi,
BGN, SRH, and bias contracts produced:

| Point | Variant | Id error | Qn p95 | n p95 |
|---|---|---:|---:|---:|
| Vg=1 V, Vd=2 V | current | 7.725% | 138.4 mV | 0.852 dex |
| Vg=1 V, Vd=2 V | Sentaurus-box + insulator | 0.981% | 9.76 mV | 0.126 dex |
| Vg=-1 V, Vd=1.1 V | current | 41.153% | 141.0 mV | 0.906 dex |
| Vg=-1 V, Vd=1.1 V | Sentaurus-box + insulator | 44.931% | 9.98 mV | 0.114 dex |

## Root-cause decision

1. **DG Id-Vd systematic high bias is primarily a DG configuration
   regression.** Restoring the earlier quantum contract reduces the 2 V Id
   error from 7.725% to 0.981% and removes most Qn/n state error.
2. **DG Id-Vg transition error should be rerun with the restored contract.**
   Existing transition spatial data and the A/B result both point to the
   self-consistent Qn branch rather than the density formula.
3. **Deep-off terminal error is a separate numerical-resolution problem.**
   The restored quantum contract greatly improves Qn/n but does not improve
   the -1 V terminal current. Exact-density substitution also has negligible
   effect. Contact-current cancellation and nonlinear/linear residual floors
   remain the leading causes.
4. **Local surface mobility and SG differences are secondary follow-ups.**
   They should be reassessed after the corrected quantum contract is used for
   the full curves, because the present self-consistent state is already on
   the wrong Qn branch.

## Next execution order

1. Promote the restored quantum contract into a separate corrected DG
   regression configuration and rerun the 21-point Id-Vd curve.
2. Rerun the 21-point DG Id-Vg curve with the same contract and re-evaluate the
   transition region separately from the three deep-off points.
3. Keep deep-off acceptance based on log-current error plus KCL resolution;
   continue contact-edge cancellation and solver precision work independently.
4. Re-run the current-carrying edge mobility/SG audit on the corrected states
   before changing mobility coefficients or Enormal discretization.

## Artifacts

- Formula replay report: `docs/validation/transportmodels_sentaurus_formula_replay_2026-08-23.json`
- Remaining-state comparison: `docs/validation/transportmodels_remaining_spatial_state_compare_2026-08-23.json`
- Quantum-contract A/B: `docs/validation/transportmodels_quantum_contract_ab_2026-08-23.json`
- Raw replay directory:
  `build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_sentaurus_formula_replay_20260823`
