# PN2D Minimal6 Sentaurus box-current staged replacement

Date: 2026-07-24

## Question

On the terminal- and KCL-closed Sentaurus box-current reconstruction, replace
the Vela inputs in the ordered chain

`QFP -> carrier density -> element mobility -> box geometry`

and determine which step closes the remaining electron and hole edge-current
error.

No production formula in `include/` or `src/` was modified.

## State and reference contract

- Topology: `mirror`.
- Bias: `-1 V`.
- Mesh: 6 physical nodes, 4 triangles, 9 global edges, and 12
  element-local edges.
- Reference current: the documented Sentaurus box-method reconstruction using
  `ReadCoefficient`, endpoint state, and element mobility.
- Electron convention: `n_plus_A_per_um`.
- Hole convention: `p_minus_A_per_um`.
- Error: `abs(log10(abs(candidate/reference)))`.
- Error statistics use the five nonzero reference edges per carrier. Exact
  zero edges remain typed and are excluded from dex statistics.

The preceding reference gate remains:

| Gate | Result |
|---|---:|
| maximum contact-current relative error | 1.114159e-7 |
| maximum total-current KCL relative error | 3.163255e-9 |
| zero-coefficient global edges | 1 and 6 |
| independent reference verification | passed |

This is a reconstructed Sentaurus box operator, not a native directed-edge
dataset. Its value as the replacement reference comes from the independent
terminal and KCL closure above.

## Ordered replacement definition

| Stage | QFP | Density | Mobility | Geometry |
|---|---|---|---|---|
| `vela_baseline` | Vela | Vela | inferred Vela production edge mobility | Vela cotangent coefficient |
| `sentaurus_qfp` | Sentaurus | Vela, frozen | Vela | Vela |
| `sentaurus_qfp_density` | Sentaurus | Sentaurus node density | Vela | Vela |
| `sentaurus_qfp_density_element_mobility` | Sentaurus | Sentaurus | Sentaurus element mobility | Vela |
| `sentaurus_qfp_density_element_mobility_geometry` | Sentaurus | Sentaurus | Sentaurus | Sentaurus `ReadCoefficient` |

The inferred Vela mobility is only the positive scalar factor of the existing
production variable-`ni` SG edge flux. The staged baseline reproduces the
integrated production edge current with a maximum relative difference of
`2.501917e-15`; it is therefore not a second current formula.

A control branch replaces the direct Sentaurus density by density recomputed
from Sentaurus `psi/QFP` and the effective intrinsic density recovered from the
Vela BGN state.

## Pooled result

| Stage | Electron median / max (dex) | Electron sign | Hole median / max (dex) | Hole sign |
|---|---:|---:|---:|---:|
| Vela baseline | 5.156389 / 5.164334 | 0.8 | 4.954712 / 5.237032 | 0.8 |
| Sentaurus QFP | 4.975977 / 5.492370 | 1.0 | 4.781812 / 5.138980 | 1.0 |
| + Sentaurus density | 0.093137 / 0.393507 | 1.0 | 0.039681 / 0.168344 | 1.0 |
| + Sentaurus element mobility | 7.71e-16 / 7.23e-15 | 1.0 | 1.54e-15 / 7.14e-15 | 1.0 |
| + Sentaurus geometry | 7.71e-16 / 7.23e-15 | 1.0 | 1.54e-15 / 7.14e-15 | 1.0 |
| recomputed-density control | 0.093141 / 0.393511 | 1.0 | 0.039685 / 0.168348 | 1.0 |

The final 10 nonzero carrier-edge values close with maximum relative error
`1.666637e-14`. The independent formula replay differs from the recorded final
stage by at most `2.635635e-15` relative.

## Paired incremental contribution

The following values are the median, over the same five edges, of
`previous error - current error`. Positive values reduce error. They are
paired statistics and need not equal differences of pooled medians.

| Replacement step | Electron (dex) | Hole (dex) |
|---|---:|---:|
| Vela baseline -> Sentaurus QFP | -0.124886 | +0.098052 |
| QFP -> QFP + density | +4.904236 | +4.769842 |
| density -> element mobility | +0.093137 | +0.039681 |
| mobility -> geometry | 0.000000 | 0.000000 |

QFP replacement alone fixes all nonzero current signs, but it does not close
the magnitude. The decisive step is using carrier densities consistent with
the replaced QFP state. Element mobility then removes the remaining
sub-0.4-dex edge variation. Geometry contributes exactly zero on this mesh.

## Per-edge error lattice

| Carrier | Sentaurus edge | Nodes | Vela baseline | QFP | QFP + density | Full closure |
|---|---:|---|---:|---:|---:|---:|
| electron | 2 | 4-5 | 0.093143 | 0.218738 | 0.093137 | 6.75e-16 |
| electron | 3 | 0-1 | 0.128067 | 0.252953 | 0.128061 | 7.71e-16 |
| electron | 4 | 1-5 | 5.164334 | 5.492370 | 0.393507 | 7.71e-16 |
| electron | 5 | 1-2 | 5.158323 | 4.979397 | 0.075161 | 7.14e-15 |
| electron | 7 | 3-5 | 5.156389 | 4.975977 | 0.039937 | 7.23e-15 |
| hole | 2 | 4-5 | 4.968081 | 4.793416 | 0.023573 | 7.14e-15 |
| hole | 3 | 0-1 | 4.954712 | 4.781812 | 0.003821 | 7.14e-15 |
| hole | 4 | 1-5 | 5.237032 | 5.138980 | 0.168344 | 1.54e-15 |
| hole | 5 | 1-2 | 0.066769 | 0.188828 | 0.066762 | 7.71e-16 |
| hole | 7 | 3-5 | 0.039687 | 0.162622 | 0.039681 | 8.68e-16 |

The central edge 1-5 is the largest remaining mobility-sensitive edge:
Vela/Sentaurus mobility differs by `0.393507 dex` for electrons and
`0.168344 dex` for holes. Across the five active edges, the mobility step's
paired median is only `0.093137 dex` and `0.039681 dex`, respectively.

## Density self-consistency control

Using the Vela BGN effective intrinsic density with the Sentaurus
electrostatic and quasi-Fermi potentials reproduces the exported Sentaurus
node densities to a maximum `4.426180e-6 dex` over all six nodes and both
carriers.

The current error with these recomputed densities differs from the direct
Sentaurus-density branch by only about `4.4e-6 dex`. Therefore the decisive
density step is not an arbitrary data substitution: it is reproduced by the
current Vela BGN equilibrium transformation once the Sentaurus potentials are
used.

## Geometry check

All 12 element-local Vela cotangent coefficients are exactly equal, in the
recorded floating-point values, to the Sentaurus `ReadCoefficient` values.
The two diagonal local supports have coefficient zero; the other coefficients
are `0.25`, `1`, or, after summing two adjacent elements on the central edge,
`2`. The maximum element-local absolute difference is `0.0`.

This directly confirms that the right-triangle/zero-dual diagonal geometry is
not the source of the remaining current discrepancy.

## Scientific conclusion

For the mirror `-1 V` state and the closed Sentaurus box reference:

1. The Vela production edge operator has been replayed exactly, so the staged
   result is not caused by a baseline implementation mismatch.
2. QFP replacement establishes the correct current direction but leaves
   approximately five decades of magnitude error on the pooled support.
3. Recomputing or replacing the carrier densities after the QFP replacement
   removes approximately `4.8-4.9 dex` and is the dominant closure step.
4. Sentaurus element mobility removes the remaining `0.04-0.09 dex` median
   error, with the central electron edge reaching `0.3935 dex`.
5. Vela and Sentaurus box geometry are identical here and contribute exactly
   zero error.
6. The full substitution closes the reconstructed Sentaurus edge current to
   floating-point precision. This supports the documented box-current algebra
   and rejects geometry or carrier-sign convention as the residual cause.

The result does not by itself authorize a production formula change. It is a
single-state operator audit, and the reference remains a terminal-closed
reconstruction rather than a native Sentaurus directed-edge observation. The
next extension, if needed, is to repeat the same sealed replay over the exact
40-state lattice.

## Evidence

- Evidence root:
  `build-release/pn2d-minimal6-sentaurus-box-staged-replacement-20260724-a`
- Stage samples: `stage_edge_samples.csv`
- Stage summary: `stage_summary.csv`
- Paired contributions: `paired_contributions.csv`
- Density control: `density_recompute_control.csv`
- Mobility comparison: `mobility_comparison.csv`
- Geometry comparison: `geometry_comparison.csv`
- Baseline cross-check: `baseline_operator_crosscheck.csv`
- Independent verification: `independent_verification.json`
- Sealed raw inputs and scripts: `raw/`
- SHA-256 ledger: `evidence_manifest.json`

The independent verifier checked all 108 stage samples, recomputed the final
Sentaurus element formula without importing the main analyzer, regenerated
all pooled medians and paired contributions, and reported zero failures.
