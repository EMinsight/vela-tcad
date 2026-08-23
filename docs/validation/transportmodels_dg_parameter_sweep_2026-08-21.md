# TransportModels DG parameter and unit audit

Work point: Vg = 1.0 V, Vd = 2.0 V. Fixed hybrid state, p1_direct Eq. 231.

Status: **pass**

## Main finding

The best tested one-factor setting is `conduction_band_narrowing_fraction=1`, with a normalized
free-node residual L1 of `0.567074`. The independently
recovered coefficient mass `1.090650673229639` gives
`1.017619` versus the explicit baseline
`1.000000`.

This scan is diagnostic, not a fitted production calibration. A lower fixed-
state residual identifies a parameter direction worth testing self-consistently;
it does not by itself prove improved terminal-current agreement.

The TDR-derived material contract reduces the residual L1 to
`0.583322` (a reduction of
`41.67%`), close to the unphysical
`Bgn2Chi=1` control at `0.567074`. This shows that
the apparent BGN-share preference is largely a proxy for base-affinity
mismatch: Si/PolySi require +22.740 mV and SiO2 requires -50.000 mV relative
to the original Vela material file. The semantic BGN fraction remains 0.5.

## Parameter mapping

| Sentaurus quantity | Vela field | Value | Unit | Mapping status |
|---|---|---:|---|---|
| QuantumPotentialParameters.gamma[electron] | `solver.electron_quantum_potential.gamma` | 3.6 | dimensionless | exact |
| QuantumPotentialParameters.theta[electron] | `solver.electron_quantum_potential.theta` | 0.5 | dimensionless | exact |
| QuantumPotentialParameters.xi/eta/nu[electron] | `fixed Eq. 231 semantics` | 1 / 1 / 0 | dimensionless | neutral defaults represented |
| eDOSMass Formula 1 at 300 K | `effective_mass_ratio` | 1.0618016171622988 | m*/m0 | frozen at 300 K |
| Formula-0 quantum coefficient mass | `coefficient_mass_ratio` | 1.0906506732296395 | m*/m0 | recovered oracle |
| Bandgap.Bgn2Chi | `conduction_band_narrowing_fraction` | 0.5 | dimensionless | Eq. 231 drive mapping |

The TransportModels `pp13_des.par` file overrides SRH lifetime parameters but
does not override the quantum-potential section, so the Silicon defaults own
the DG parameter values. The two electron mass roles are intentionally kept
separate: DOS mass changes the material drive, while coefficient mass changes
the gradient coefficient.

## One-factor scan

| Variant | Residual L1 / baseline | Max residual / baseline | Substrate L1 share | Max node |
|---|---:|---:|---:|---:|
| conduction_band_narrowing_fraction=1 | 0.567074 | 0.997692 | 43.91% | 2 |
| corrected_material_contract | 0.583322 | 1.017215 | 46.04% | 2 |
| gamma=4 | 0.935223 | 1.000000 | 66.90% | 2 |
| coefficient_mass_ratio=1 | 0.962283 | 1.000000 | 67.83% | 2 |
| effective_mass_ratio=1.09065067 | 0.979714 | 1.000000 | 68.40% | 2 |
| baseline | 1.000000 | 1.000000 | 69.04% | 2 |
| coefficient_mass_ratio=1.09065067 | 1.017619 | 1.000000 | 69.58% | 2 |
| theta=0.25 | 1.019100 | 1.084972 | 67.75% | 2 |
| theta=0.75 | 1.026704 | 0.915028 | 67.63% | 2 |
| effective_mass_ratio=1 | 1.045448 | 1.000000 | 70.39% | 2 |
| theta=0.1 | 1.054908 | 1.135955 | 65.62% | 2 |
| gamma=3.2 | 1.081140 | 1.000000 | 71.36% | 2 |
| coefficient_mass_ratio=1.2 | 1.084492 | 1.000000 | 71.45% | 2 |
| conduction_band_narrowing_fraction=0 | 1.525269 | 1.002292 | 78.71% | 2 |

All variants use identical mesh, state, contacts, discretization, temperature,
and carrier models. Only the named DG parameter changes.
