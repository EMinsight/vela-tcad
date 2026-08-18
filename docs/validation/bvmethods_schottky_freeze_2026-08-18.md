# BVmethods non-transient and simple Schottky freeze

## Result

The simple Schottky workflow is closed at **PASS**. The Sentaurus BVmethods
non-transient reference set is frozen at **PASS_WITH_CONTINUATION_PENDING**:
four Vela mappings are accepted, while the NMOS pseudo-arclength mapping is
explicitly not claimed complete. Sentaurus transient node 8 is outside scope.

| Scope | Sentaurus reference | Vela result | Status |
| --- | ---: | ---: | --- |
| BVmethods ABA Poisson | 5.305525633 V | operator reference | frozen |
| BVmethods coupled ABA/IIC | 6.377494278 V | path IIC pass; current IIC within 3% | pass |
| BVmethods external resistor | 6.379791636 V | 6.401809065 V | pass, 0.345% |
| BVmethods voltage-to-current | 6.383184201 V | 6.401911823 V | pass, 0.293% |
| BVmethods Continuation | 6.383727169 V | no usable NMOS arc curve | pending |
| Simple Schottky I-V | 0.01--1 V, 24 points | max 0.478652 dex | pass |

## BVmethods freeze boundary

The five compact Sentaurus curves are generated from archived O-2018.06-SP2
PLT files and sealed by SHA-256 in
`reference_tcad/bvmethods_sentaurus2018/bvmethods_nontransient_validation_20260817.json`.
The transient deck is retained as source material but is not run or accepted.

The NMOS Continuation investigation did not change the mesh, contacts,
materials, mobility, SRH, impact-ionization model, or coefficients. Both the
full-physics deck and a no-Enormal control accepted only a roughly `1.5e-10 V`
advance before reaching the minimum arclength step. A second trial used two
already accepted full-physics states at 6.000000 V and 6.056459 V to initialize
the first tangent by a secant; after eight minutes it still had no accepted
point beyond the anchor. The evidence therefore points to state scaling and
bordered-corrector robustness, not to a missing physical model. Details are in
`reference_tcad/bvmethods_sentaurus2018/continuation_diagnostic_20260818.json`.

## Simple Schottky closeout

`scripts/run_schottky_reference_workflow.py` creates a clean stage-A voltage
sweep to 0.82 V, restarts that exact state into stage-B pseudo-arclength, and
merges the curve through 1 V. It gates all 24 nonzero Sentaurus points and
checks terminal KCL near 0, 0.4, and 1 V. The maximum current-shape error is
0.478652 dex; at 1 V the Vela current is `1.08262e-4 A/um` versus
`1.18824e-4 A/um` in Sentaurus, an 8.89% difference. No image-force lowering,
tunnelling, series resistance, AC, or high-field model is enabled.

## Reproduction and acceptance artifacts

- `scripts/freeze_bvmethods_nontransient.py`: rebuilds the five BVmethods
  reference curves and acceptance ledger.
- `scripts/prepare_bvmethods_nmos_continuation.py`: prepares the physics-frozen
  NMOS continuation diagnostic deck.
- `scripts/run_schottky_reference_workflow.py`: runs and audits the complete
  two-stage Schottky comparison.
- `tests/regression/test_bv_schottky_frozen_ledgers.py`: prevents either
  checked-in ledger from silently losing its accepted scope or status.

Generated state files, nonlinear histories, TDR/PLT files, and full simulation
directories remain ignored build artifacts.
