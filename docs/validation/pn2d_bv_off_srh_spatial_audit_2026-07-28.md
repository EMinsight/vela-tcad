# PN2D avalanche-off SRH spatial baseline

Date: 2026-07-28

## Frozen configuration

- Starting commit: `fa1c343`
- BV mobility doping basis: `net_doping`
- Forward-IV mobility doping basis: `cell_reconstructed_total_impurity`
- Impact ionization: disabled
- SRH and Old Slotboom BGN: enabled
- Bias points: `0, -1, ..., -20 V`
- Sentaurus source baseline:
  `reference_tcad/pn2d_sentaurus2018_coarse7x3/source`
- Physical mesh hash in both simulators:
  `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e`
- Physical mesh, doping, contacts, temperature, device depth, and Newton
  residual tolerances: unchanged
- The independent global continuity-closure acceptance gate was tightened to
  `1e-5` so Task 1 cannot accept a residual-closed but source/contact-open
  point.

## Task 1 acceptance evidence

The aligned artifacts are under:

`build-release/pn2d-bv-off-srh-spatial-audit-coarse-baseline-20260728/`

| Gate | Result | Limit | Status |
| --- | ---: | ---: | --- |
| Vela convergence | 21/21 | 21/21 | pass |
| Sentaurus convergence/TDR coverage | 21/21 | 21/21 | pass |
| max electron source/contact closure | `7.482e-7` | `1e-5` | pass |
| max hole source/contact closure | `3.561e-9` | `1e-5` | pass |
| max total terminal closure | `3.762e-23 A/um` | `1e-20 A/um` | pass |
| exported-source reintegration error | `4.893e-15` | `1e-6` | pass |
| constant/linear triangle integration tests | 3/3 | `1e-12` relative | pass |

The VTK diagnostic path now writes floating-point values at
`max_digits10`. Before that repair, six-digit VTK rounding alone produced a
`1.25e-6` source reintegration error and triggered the Task 1 stop condition.

## Terminal-current anchors

| Bias (V) | Vela (A/um) | Sentaurus (A/um) | `log10(abs(Vela/Sentaurus))` |
| ---: | ---: | ---: | ---: |
| -1 | `3.098001e-17` | `3.098171e-17` | `-2.385e-5` |
| -5 | `3.101277e-17` | `3.101630e-17` | `-4.939e-5` |
| -10 | `3.106216e-17` | `3.106658e-17` | `-6.170e-5` |
| -15 | `3.113177e-17` | `3.113760e-17` | `-8.124e-5` |
| -20 | `5.143587e-17` | `5.142938e-17` | `+5.477e-5` |

These values replace the former non-paired reference. They show that the
coarse7x3 Sentaurus source and Vela mesh are already paired to much better
than `0.05 dex`. The former shape discrepancy must not be used as evidence
for the new baseline.

## Spatial data contract

The report exports native Vela and Sentaurus node and triangle tables,
cumulative absolute-source profiles, source centroids and 10/50/90%
positions, positive generation and negative recombination integrals, and the
fixed depletion indicator

`(n + p) < 0.1 * max(abs(net_doping), effective_intrinsic_density)`.

The junction-normal coordinate is `x - x_junction`, where `x_junction` is
inferred from the signed doping zero crossing. Sentaurus does not provide a
native effective-intrinsic-density field in these TDR files; this is recorded
as unavailable and is not reconstructed or silently zero-filled. Native
Sentaurus `psi`, `n`, `p`, doping, and `SRHRecombination` are available at all
21 points.

Every figure is mapped to its source CSV in `report/report_manifest.json`.
The run manifest records the Git dirty flag, configuration and input hashes,
TDR hashes, model switches, and the frozen BV/IV mobility-basis choices.

## Task 1 decision

Task 1 passes. The exported Vela source exactly reproduces the solver source
well inside the stop-condition threshold, so cross-simulator same-state
decomposition may proceed. No scientific root-cause classification is made
at this gate.
