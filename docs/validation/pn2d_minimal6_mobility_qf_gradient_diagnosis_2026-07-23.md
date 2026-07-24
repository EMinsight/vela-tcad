# PN2D Minimal6 mobility and quasi-Fermi-gradient diagnosis

> Quantitative correction: the direction conclusion is unchanged, but the
> original numeric tables used pre-fix Vela restart-unit inputs. Corrected
> values and evidence roots are in
> `pn2d_minimal6_state_unit_and_transport_reaudit_2026-07-23.md`.

Date: 2026-07-23

Scope: 40 exact Sentaurus/Vela states, two topologies, two carriers

Production-formula changes: none

## Conclusion

The earlier approximately 90 degree aggregate direction error is a carrier-sign pooling artifact, not a mobility-induced rotation. On the exported Sentaurus quasi-Fermi-potential convention, both carrier-current vectors satisfy the data-supported local orientation

\[
\mathbf{J}_n \parallel -\nabla \phi_{Fn}, \qquad
\mathbf{J}_p \parallel -\nabla \phi_{Fp}.
\]

Applying the negative sign to both exported carrier potentials reduces the pooled median angle from 89.965351 degrees to 0.010067 degrees. Identity is the best coordinate transform for electrons and holes; axis swaps and 90 degree rotations give 90 degrees, while full negation gives approximately 180 degrees.

Mobility affects the inverted gradient magnitude but cannot explain the old direction error. The decisive magnitude effect is projection order: locally forming `J/(q n mu)` and then projecting is materially better than separately projecting `J`, `n`, and `mu` before division.

## Edge-mobility comparison

All mobilities are in m2/(V s). Errors are absolute log10 errors against the exported Sentaurus mobility on the same endpoint-mean edge support.

| Carrier | Branch | N | Median mobility | Median abs error (dex) | P95 abs error (dex) | Median relative error |
|---|---|---:|---:|---:|---:|---:|
| electron | Sentaurus exported | 360 | 0.0153903 | 0 | 0 | 0 |
| electron | Vela Masetti, Sentaurus state | 360 | 0.0294850 | 0.102586 | 0.983253 | 0.266445 |
| electron | Constant 1417 cm2/(V s) | 360 | 0.1417000 | 0.964305 | 1.171215 | 8.214804 |
| electron | Vela Masetti, native state | 360 | 0.0299165 | 0.100664 | 0.983338 | 0.260852 |
| hole | Sentaurus exported | 360 | 0.0110031 | 0 | 0 | 0 |
| hole | Vela Masetti, Sentaurus state | 360 | 0.0189197 | 0.086174 | 0.664796 | 0.219477 |
| hole | Constant 470.5 cm2/(V s) | 360 | 0.0470500 | 0.631178 | 0.819144 | 3.278709 |
| hole | Vela Masetti, native state | 360 | 0.0191653 | 0.084403 | 0.664794 | 0.214515 |

The same-state Masetti branch is much closer to Sentaurus than the constant control, but its tails remain large, especially the electron P95 near 0.98 dex. A single global mobility scale is therefore not supported.

## Quasi-Fermi-gradient inversion

The table reports the median absolute log10 magnitude error. The direction angles shown for cell-local inversion use the corrected negative sign for both exported carrier potentials.

| Support and order | Mobility | Electron median / P95 (dex) | Hole median / P95 (dex) | Electron / hole median angle |
|---|---|---:|---:|---:|
| edge, ratio after projection | Sentaurus exported | 0.210796 / 12.494413 | 0.230585 / 12.387723 | tangent sign agreement 1.0 / 1.0 |
| edge, local inversion then projection | Sentaurus exported | 0.128209 / 0.504766 | 0.127876 / 0.517886 | tangent sign agreement 1.0 / 1.0 |
| cell, ratio after projection | Sentaurus exported | 5.779672 / 12.572908 | 5.687322 / 12.461543 | 0.006020 / 0.006184 deg |
| cell, local inversion then projection | Sentaurus exported | 0.142764 / 0.675465 | 0.149314 / 0.688801 | 0.009142 / 0.011362 deg |
| cell, local inversion then projection | Vela Masetti, Sentaurus state | 0.672560 / 1.217038 | 0.487891 / 1.002407 | 0.009802 / 0.011902 deg |
| cell, local inversion then projection | constant control | 1.197477 / 1.783458 | 0.871270 / 1.449292 | 0.009133 / 0.011352 deg |

Changing the mobility branch changes the cell-local direction by at most approximately 0.00067 degrees. In contrast, forming a ratio after projection adds more than 5 dex of cell magnitude error. Averaging and division do not commute across this strongly nonuniform junction.

Eighty of 360 carrier-edge rows have zero reference tangent quasi-Fermi gradient, leaving 280 magnitude-valid rows per carrier. They remain counted for support auditing but are excluded from log-magnitude summaries.

## Orientation controls

| Carrier | Identity median | Negate median | Swap XY median | Rotate CW/CCW median |
|---|---:|---:|---:|---:|
| electron | 0.009142 deg | 179.990858 deg | 90.000000 deg | 90.000000 deg |
| hole | 0.011362 deg | 179.988638 deg | 90.000000 deg | 90.000000 deg |

The legacy mixed-sign pooling reproduces 89.965351 degrees. The corrected same-sign pooling gives 0.010067 degrees. This rejects mobility, XY swapping, and coordinate rotation as explanations of the former near-orthogonal result.

## Method and limitations

The diagnostic reproduces the production Masetti low-field formula and quasi-Fermi-gradient field limiter in SI units, while leaving `include/` and `src/` unchanged. Sentaurus node mobility is mapped to edges by endpoint means. Cell gradients use the affine P1 triangle gradient. Two inversion orders are kept separate:

1. project current, density, and mobility independently, then divide;
2. form the local node quantity `-J/(q n mu)`, then project it.

Sentaurus current in this bundle is node-exported, not an internal directed-edge flux. The local drift identity `J/(q n mu)` is not a discrete Scharfetter-Gummel inverse, so the remaining approximately 0.13--0.15 dex median magnitude gap is not sufficient evidence for a production mobility or current-formula change.

## Recommended next discriminating test

Export or reconstruct Sentaurus current/flux on exactly the same directed edges, then invert the discrete Scharfetter-Gummel relation using endpoint densities and quasi-Fermi-potential differences. This is the next test capable of separating current semantics, support mapping, and mobility magnitude without fitting an empirical scale.

## Reproducible artifacts

- Diagnostic CLI: `scripts/diagnose_pn2d_minimal6_mobility.py`
- Implementation: `scripts/pn2d_minimal6_diagnostics/mobility_diagnosis.py`
- Regression tests: `tests/regression/test_pn2d_minimal6_mobility_diagnosis.py`
- Authoritative output: `build-release/pn2d-minimal6-mobility-diagnosis-20260723-b`
- Independent deterministic repeat: `build-release/pn2d-minimal6-mobility-diagnosis-20260723-c`
- Portable report: `build-release/pn2d-minimal6-mobility-diagnosis-20260723-b/report.html`
