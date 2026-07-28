# PN2D avalanche-off continuity closure validation (2026-07-28)

## Scope

This validation covers the source-aware continuity convergence repair for the
coarse7x3 PN2D reverse-bias case with impact ionization disabled.

## Root cause and repair

The TCAD-unit implementation of
`PhysicalUnitSystem::continuitySourceIntegralFactor()` used
`area / mobility = 1e-8`.  That expression did not convert a volumetric rate
to the same per-device-depth units as an SG edge-flux line integral.  The
correct conversion is

```text
concentration_scale * area_scale
--------------------------------- = 1e-4
current_density_scale * length_scale
```

The old factor understated every volumetric continuity source by `1e4`.
The repair also:

- enforces per-carrier global source/contact-flux closure;
- scales source-dominated continuity rows independently;
- solves high-reverse-bias quasi-Fermi potentials as majority-contact
  referenced increments;
- preserves those increments through terminal-current and VTK diagnostics.

## Avalanche-off results

Artifacts:

- `build-release/pn2d-bv-off-rootcause-20260728/source_scale_fix_full20/`

| Bias (V) | Vela off (A/um) | Sentaurus off (A/um) | Vela/Sentaurus |
|---:|---:|---:|---:|
| -1  | 3.0980e-17 | 8.70e-18 | 3.56 |
| -5  | 3.1013e-17 | 2.73e-17 | 1.14 |
| -10 | 3.1062e-17 | 4.21e-17 | 0.738 |
| -15 | 3.1132e-17 | 5.26e-17 | 0.592 |
| -20 | 5.1436e-17 | 6.16e-17 | 0.835 |

At each nonzero anchor, independently integrated SRH generation agrees with
the summed electron and hole contact fluxes to approximately `2e-6` relative
or better.  The sum of both terminal total currents closes to approximately
`4e-23 A/um` or better.

The original two-order-of-magnitude leakage deficit is therefore removed.
The remaining shape difference from -1 V through -15 V is physical/model or
spatial-discretization work, not a continuity-source loss or convergence
false positive.

## Forward-IV regression

The 0-20 V `cell_reconstructed_total_impurity` candidate converged at all 201
bias points.  Relative to the pre-repair curve:

| Bias (V) | Relative current change |
|---:|---:|
| 0.5 | +0.2593% |
| 1 | -0.00179% |
| 5 | -0.000161% |
| 10 | -0.0000724% |
| 15 | -0.0000463% |
| 20 | -0.0000338% |

The source-scale repair does not materially disturb the established forward
production curve above the low-current transition region.
