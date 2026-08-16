# SingleDevice Eq. 231 region-side trace implementation (2026-08-14)

> 2026-08-15 correction: the direct SingleDevice box file represents signed
> `AverageBoxMethod` circumcentric measures, not the positive mixed-Voronoi
> obtuse fallback.  The corrected audit and endpoint result are documented in
> `singledevice_eq231_averagebox_polyoxide_central_closure_2026-08-15.md`.

## Outcome

The exported `eQuantumPotential` value is now treated as an owner/output
trace.  The experimental `sentaurus_box` assembly reconstructs independent
Si, PolySilicon, and insulator-side fitted traces at a shared vertex and uses
material-side reaction weights with a consistent diagonal Jacobian.

This closes the previously identified interface implementation gap, but it
does not complete the self-consistent endpoint gate.  After the interface
rows were reduced, the leading fixed-state discrepancy moved to ordinary
Silicon bulk rows.  Both imported endpoints still move by about `1.000394 V`,
so the 21-point curves remain gated.

## Direct evidence

- VM `MeasureCoefficients.debug` was copied read-only and compared with the
  local measures.  At the probed SingleDevice interface cells the ratio is
  one to floating-point precision, excluding a box-volume mismatch.
- Two independent oxide/PolySilicon Jacobian probes and 71 adjacent oxide
  rows give a common insulator half-jump offset near `0.02012`.
- Adjacent pure-material rows recover a Silicon-side offset of about
  `-4.60e-5` and a PolySilicon-side offset of about `0.0026675`.
- A 65-node PolySilicon/SiO2 fit requires separate material-side reaction
  traces; a single global multiplier is rejected.

The supplied GSS book supports the bulk fitted operator but does not disclose
the missing Sentaurus interface algebra.  Its Eq. 9.128 is the same
exponential/second-order DG flux used here and explicitly notes its local
nonconservation.  Section 9.11.9 says that an explicit SiO2/multimaterial
global DG solve keeps the exported quantum potential continuous.  This is
consistent with treating the TDR value as a continuous owner/output trace,
while reconstructing internal material-side numerical traces during assembly.

## Implementation

- Replaced the invalid global interface reaction multiplier by four explicit
  side weights: Silicon, PolySilicon, insulator-at-Silicon, and
  insulator-at-PolySilicon.
- Added three explicit fitted half-jump offsets: insulator, Silicon, and
  PolySilicon.
- Applied the same side factor to residual and analytic Jacobian.
- Corrected all-material gate-contact quantum handling: an Ohmic gate contact
  is Dirichlet; only the semiconductor-only truncation path skips a metal-gate
  contact.
- Added fixed-state fit/audit scripts and a manufactured two-material restart
  test.  All controls are opt-in under `sentaurus_box` and default to neutral
  values.

## Verification

- Density-gradient tests: 86 assertions in 21 cases pass.
- Newton density-gradient parser tests pass; the side controls are parsed and
  range-checked.
- Linear fixed state: leading interface residual is about `0.643`; the global
  maximum is `2.470`, now at ordinary Silicon node 2101.
- Saturation fixed state: global maximum is `2.497`, also at Silicon node 2101.
- Linear endpoint: inner quantum solve converges in 22 iterations, but raw
  imported-state change is `1.000394 V`; endpoint fails.
- Saturation endpoint: inner quantum solve converges in 22 iterations, but raw
  imported-state change is `1.000394 V`; endpoint fails.

## Decision

The region-side interface-trace task is implemented and isolated.  The next
root-cause target is no longer the interface, gate boundary, box measure, or
ordinary damping.  It is the Formula-0 fitted bulk-Silicon source/reaction
closure, beginning with nodes 2101, 2100, and 2044 under both endpoint states.
Do not run the complete 21-point Id-Vg curves until that fixed-state bulk gate
and both endpoint gates pass.
