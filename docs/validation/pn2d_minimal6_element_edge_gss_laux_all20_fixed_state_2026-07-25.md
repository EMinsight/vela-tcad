# PN2D Minimal6 all-20 element-edge fixed-state audit

Date: 2026-07-25

## Outcome

Task 8 passed for the exact 40-state matrix:

- topologies: `mirror`, `sketch`
- biases: -1 V through -20 V
- state variables imported from Sentaurus: `psi`, `phin`, `phip`
- quantities recomputed by Vela: `n`, `p`, mobility, local SG edge
  currents, GSS/Laux cell vectors, Van Overstraeten coefficients, and
  element-vertex avalanche source integrals

The opt-in element-edge diagnostic closes the documented Sentaurus box
operator over the complete bias range. This is evidence for the diagnostic
operator and does not change the production default.

## Sentaurus regeneration

The remote Sentaurus O-2018.06-SP2 run completed both topologies. Each
topology contains exactly 20 observed targets and no missing or extra bias.
The raw manifest is:

`build-release/pn2d-minimal6-element-avalanche-all20-runtime-20260725/manifest.json`

For each topology, the runtime log contains:

| Record | Per state | Per topology | Both topologies |
| --- | ---: | ---: | ---: |
| physical plus electrode vertices | 10 | 200 | 400 |
| elements | 4 | 80 | 160 |
| element-vertex measures | 12 | 240 | 480 |
| local directed element edges | 12 | 240 | 480 |
| device integrals | 1 | 20 | 40 |

## Documented Sentaurus operator replay

The deterministic analysis root is:

`build-release/pn2d-minimal6-element-avalanche-all20-runtime-20260725/analysis`

Across 80 topology-bias-carrier samples, the device source-integral errors
against native Sentaurus avalanche integrals are:

| Reconstruction | Median error (dex) | Maximum error (dex) |
| --- | ---: | ---: |
| GSS/Laux edge-volume weighted | 0.0000473942 | 0.000429237 |
| active-edge exact control | 0.0000473942 | 0.000429237 |
| Charon Whitney cell average | 0.00121779 | 0.00285592 |
| Genius tangent least squares | 0.00149059 | 0.00365110 |

The GSS/Laux and active-edge controls are identical on this four-right-
triangle mesh because the diagonal box partial areas are zero.

## Vela fixed-state replay

The deterministic fixed-state root is:

`build-release/pn2d-minimal6-element-edge-gss-laux-fixed-state-all20-20260725`

The independent verifier checked 480 element-vertex rows, 40 states, 160
zero diagonal partial volumes, and the exact manifest state lattice.

| Quantity | Electron median (dex) | Electron max (dex) | Hole median (dex) | Hole max (dex) |
| --- | ---: | ---: | ---: | ---: |
| directed SG edge current | 0.0363836 | 1.06335 | 0.0338998 | 0.710160 |
| GSS/Laux cell current vector | 0.0337634 | 0.220792 | 0.0318710 | 0.0847706 |
| electric-field avalanche coefficient | 1.02e-12 | 1.55e-9 | 1.68e-12 | 2.56e-9 |
| accumulated node source | 0.0406078 | 0.232730 | 0.0379011 | 0.0902511 |
| device source integral | 0.0548177 | 0.197052 | 0.0517363 | 0.0767377 |

The element source identity,

`qG = alpha * abs(J) * element_vertex_measure`,

closes to a maximum relative error of `4.10e-16`.

The large maximum directed-edge errors occur on tiny-current central edges.
They do not dominate the reconstructed cell vector or the integrated
avalanche source. The all-state medians improve relative to the sealed
-1/-10/-20 V subset.

## Scientific decision

The complete fixed-state evidence supports all of the following:

1. Sentaurus Van Overstraeten coefficients in this experiment are driven by
   the element electric-field magnitude, while imported QFP values determine
   carrier density and SG current.
2. On the Minimal6 right triangles, a zero diagonal box partial area makes
   that edge numerically inactive in the GSS/Laux vector; this is why the
   GSS/Laux and two-active-edge controls are exactly identical. The general
   implementation still supports all three edges for non-right triangles.
3. The dominant fixed-state avalanche-source discrepancy was the triangle
   current support and vector reconstruction, not the Van Overstraeten
   coefficient formula, mobility formula, or QFP sign convention. This does
   not classify the remaining self-consistent QFP/current mismatch.

Task 9 must now test the same operator in a self-consistent 40-state sweep.
No production default should change until that nonlinear sweep and the
general-mesh regressions pass.
