# PN2D Minimal6 Task 8 corrected validation and physics error audit (2026-07-15)

## Conclusion

Tasks 1-8 are complete for the Minimal6 diagnostic scope. The corrected Task 8
comparison contains 40 exact common checkpoints (sketch/mirror x -1..-20 V),
reaches -20 V for both solvers and has no rejected common checkpoint. These
curves remain diagnostic sweeps, not physical BV curves.

The dominant discrepancy is not the electrostatic potential, electric field or
the Van Overstraeten ionization coefficient. In the semiconductor interior,
Vela and Sentaurus electrostatic potential agrees to at most 2.31e-11 V and the
maximum electric-field ratio is 1.001-1.014. At -12..-20 V the peak electron
ionization coefficient ratio is 1.015-1.040. The dominant state error is the
absolute electron and hole quasi-Fermi-potential branch: the interior error is
0.30-0.37 V, producing carrier-density excesses of 4.7e4-8.0e5 through the
Boltzmann exponential and terminal-current ratios of 2.24e3-8.34e3.

The code audit identifies a physical unit inconsistency in the coupled
continuity residual as the leading implementation-level cause. SG edge
transport is assembled with micrometre geometry while volumetric recombination
and generation are multiplied by an area in um^2, without the matching
line-to-centimetre and area-to-square-centimetre factors. This overweights a
volumetric source relative to transport by 1e4. A separate 1e-8 omission in the
diagnostic source export was corrected; it affected the reported source unit,
not the nonlinear state.

## Task 1-8 completion ledger

| Task | Status | Evidence |
| --- | --- | --- |
| 1 | Complete | Six exact extended-field states and independently sealed recovery package. |
| 2 | Complete | Ledger contracts, units and deterministic report schemas. |
| 3 | Complete | Independent geometry, physics, alpha, intrinsic-density, integration and support conversions. |
| 4 | Complete | Fixed-state report validated and retained the evidence-limited `insufficient_data` conclusion. |
| 5 | Complete | Five reviewed PNG/PDF figure pairs and manifest. |
| 6 | Complete | Vela and Sentaurus drivers with retained endpoints and failures. |
| 7 | Complete | Corrected comparison regenerated from the final Vela and Sentaurus roots. |
| 8 | Complete | Corrected sweeps, comparison, physics audit, regression verification and this validation record. |

Authoritative corrected roots:

- Vela: `build-release/pn2d-minimal6-vela-task8-source-unit-corrected-20260715`
- Sentaurus corrected current manifest: `build-release/pn2d-minimal6-sentaurus-task8-current-corrected-20260715`
- Sentaurus full-field exports used for node analysis: `build-release/pn2d-minimal6-sentaurus-task8-final-r2-20260715`
- Comparison: `build-release/pn2d-minimal6-comparison-task8-source-unit-corrected-20260715`
- Self-consistent node analysis: `build-release/pn2d-minimal6-self-consistent-physics-task8-20260715`
- Fixed-state node figures: `build-release/pn2d-minimal6-node-physics-task8-20260715`
- Interactive report: `build-release/pn2d-minimal6-physics-error-report-task8-20260715/report.html`

Earlier `final-r2`, `final-r3`, `physics-corrected` and `qf-pinned-r2` roots are
intermediate or superseded and must not be cited as the Task 8 result.

## Corrections applied before the final comparison

1. Sentaurus TDR coordinates are metres; topology inputs now convert x and y to
   micrometres with a factor of 1e6 and declare `coordinate_unit: um`.
2. For this two-dimensional comparison, Sentaurus `ContactCurrentFlux` in A is
   compared numerically with the Vela A/um quantity. The earlier additional
   factor of 1e-4 was removed.
3. A reconstructed edge source has the raw dimensions
   `alpha[cm^-1] * flux[cm^-2 s^-1] * area[um^2]`. It is multiplied by 1e-8 to
   publish `s^-1/cm`.
4. Minority-electron contact relaxation is explicitly disabled in the replay
   template. This pins contact quasi-Fermi values but does not remove the
   interior quasi-Fermi or terminal-current discrepancy.

## Self-consistent comparison

The table shows the sketch topology; mirror results are numerically identical
at the quoted precision.

| Bias (V) | max interior abs psi error (V) | mean abs phin error (V) | mean abs phip error (V) | max E ratio | max alpha_n ratio | geometric n ratio | terminal current ratio | corrected source ratio |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| -1 | 2.31e-11 | 0.3032 | 0.2955 | 1.0143 | not meaningful at approximately 1e-53/cm | 6.33e4 | 3.18e3 | 1.06e4 |
| -12 | 4.23e-12 | 0.3318 | 0.3308 | 1.0020 | 1.0396 | 1.91e5 | 2.24e3 | 1.13e3 |
| -19 | 2.84e-12 | 0.3470 | 0.3468 | 1.0013 | 1.0164 | 3.45e5 | 3.66e3 | 1.80e3 |
| -20 | 1.05e-12 | 0.3686 | 0.3685 | 1.0012 | 1.0148 | 7.96e5 | 8.34e3 | 4.09e3 |

The maximum all-node potential error is 0.01295 V and occurs at contact/built-in
boundary nodes. It does not represent the semiconductor-interior agreement.

### Electrostatic potential and electric field

The interior potential match and the 0.1-1.4% peak-field difference exclude
Poisson electrostatics and field reconstruction as the dominant source of the
three-to-four-order terminal-current error.

### Quasi-Fermi potentials and carrier density

The carrier relations used by Vela are

```text
n = ni * exp((psi - phin) / Vt)
p = ni * exp((phip - psi) / Vt)
```

Therefore a quasi-Fermi offset `Delta phi` causes an exponential density factor
`exp(abs(Delta phi)/Vt)`. The observed 0.30-0.37 V offsets are sufficient to
explain the 1e5-scale carrier excess. At -19 V, for example, the mean interior
quasi-Fermi split is 0.01662 V in Vela and 0.71046 V in Sentaurus, while the
geometric carrier ratios are 3.45e5 for electrons and 3.43e5 for holes.

### Impact ionization

For the configured Van Overstraeten model, each carrier coefficient follows

```text
alpha(E) = gamma(T) * a * exp(-gamma(T) * b / abs(E))
```

with the configured low/high-field branch parameters. The production formula
and the independent reconstruction agree. Because E agrees, peak alpha also
agrees within about 1.5-4.0% at -12..-20 V. The -1 V ratio is not physically
informative because both coefficients are effectively zero.

Impact generation based on current density is

```text
G_II = (alpha_n * abs(Jn) + alpha_p * abs(Jp)) / q
```

Thus the remaining generation mismatch is inherited primarily from the current
and carrier-state mismatch, not from alpha(E).

## Continuity-residual formula audit

In `src/equation/CoupledDDAssembler.cpp`, SG edge transport uses a coefficient
proportional to

```text
mu * Vt * 1e4 * (couple / h)
```

where `h` and `couple` are expressed using micrometre mesh geometry. The
volumetric rows add `R * vol` and subtract `G * vol`, where `vol` is in um^2.
To integrate a cm-based current density and a cm^-3 volumetric rate over a
two-dimensional micrometre mesh, the geometric measures require

```text
edge length: 1 um = 1e-4 cm
cell area:   1 um^2 = 1e-8 cm^2
```

Using neither conversion leaves the source-to-transport relative scale too
large by `1e-8 / 1e-4 = 1e-4`; equivalently, the assembled source is overweighted
by 1e4 relative to transport. This is the strongest code-level explanation for
the wrong quasi-Fermi branch and 1e3-1e4 current discrepancy. Existing
assembler-equivalence tests compare two implementations that share this unit
convention, so they do not establish equivalence to a physical cm-based legacy
residual when volumetric sources are active.

## Recommended next correction and acceptance gates

Do not tune ionization coefficients or Poisson electrostatics. First add a
legacy-SI physical-equivalence test for continuity rows with nonzero
recombination and impact generation, then introduce explicit line and area
measure conversions in all coupled residual and Jacobian paths. Re-run the same
40-point matrix and require:

1. potential and electric-field agreement to remain unchanged;
2. interior phin/phip error and carrier ratios to collapse;
3. terminal current and integrated generation ratios to approach unity;
4. analytic and finite-difference Jacobian paths to agree;
5. sketch/mirror symmetry and all existing regressions to remain valid.

## Validation status and limitations

The numerical claims above were recomputed from the exact self-consistent node
exports at -1, -12, -19 and -20 V for both topologies, then reconciled against
the 40-row sweep comparison. The HTML report passed structural packaging and
data-query validation. Browser interaction could not be verified because no
Chromium executable is installed in the local environment. No core residual
formula has been changed in this task; the proposed continuity correction is a
root-cause conclusion requiring its own test-driven implementation task.

## Fresh verification

Overall assessment: **share with caveats**. The analysis is reproducible and
the principal calculations are cross-checked, but the continuity-residual cause
is a code-audit inference until the proposed physical-equivalence test and
correction are implemented.

- 10 Minimal6 Python regression modules: 147 tests in 80.560 s, PASS.
- test_fixed_state_operator_audit: 4 cases / 67 assertions, PASS.
- test_impact_ionization: 40 cases / 511 assertions, PASS.
- test_cell_reconstructed_avalanche: 14 cases / 83 assertions, PASS.
- Final comparison: 42 rows total, comprising 40 exact common -1..-20 V
  records and two explicit Sentaurus-only 0 V records.
- Self-consistent analysis: eight rows (two topologies x four biases), with
  all three input roots recorded in self_consistent_analysis.json.
- HTML presentation: title, source tooltips, terminal-current metric, tables
  and embedded data are present. Structural validation passed; browser
  interaction remains unverified because Chromium is unavailable.
- git diff --check: PASS.