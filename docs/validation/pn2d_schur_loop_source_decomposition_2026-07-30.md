# PN2D effective Schur-loop source decomposition

Date: 2026-07-30

Outcome:
`transport_and_avalanche_independently_sustain_reversal`

Parent outcome: `bidirectional_poisson_qfp_closed_loop_cause`

Production defaults changed: no

Task 8 authorized: no

## Contract

This observation-only diagnostic decomposes the effective Poisson-QFP loop

```text
K = C A^-1 B

A = J_psi_psi
B = J_psi_qfp
C = J_qfp_psi
D = J_qfp_qfp
S = D - K
```

at the exact `-19.7/-19.8 V` Vela avalanche-on states. It preserves the same
production Jacobian, scaling, constraints, physical state, and frozen
Sentaurus-QFP substitution residual used by the parent cross-block audit.

The diagnostic provides:

- four carrier paths: electron/hole residual rows from electron/hole QFP
  columns;
- three additive `C` model components:
  transport/boundary, SRH/Auger, and SG avalanche;
- the corresponding `C_component A^-1 B` matrices;
- model leave-one-out and model-only Schur counterfactuals;
- every matrix entry's row/column node, coordinates, carrier, value, and sign;
- raw and row/column-L2-equilibrated singular-value condition estimates; and
- double-symmetric directional finite differences for `B` and `C`.

The boundary portion of the transport/boundary component has zero nonzero
`C` entries on all six contact rows at both biases. Its observed cross-block
support is therefore interior transport, not a direct contact-row coupling.

## Causal result

| Bias | Full loop cosine | Transport-only cosine | Avalanche-only cosine | SRH/Auger-only cosine |
|---:|---:|---:|---:|---:|
| -19.7 V | -0.44177 | -0.15279 | -0.24273 | +0.75664 |
| -19.8 V | -0.45824 | -0.15546 | -0.24545 | +0.75108 |

Removing avalanche leaves the same adverse direction as the transport-only
case. Removing transport leaves the same adverse direction as the
avalanche-only case. Thus transport and avalanche are two independently
sufficient closed-loop paths at these states. SRH/Auger alone retains the
positive carrier-block direction and is excluded as a material cause.

This result does not mean that either production term should be removed.
Each model-only experiment changes only the diagnostic Schur loop and does
not define a physically closed solver.

## Magnitude, carrier, sign, and spatial support

| Bias | Transport loop L2 | Avalanche loop L2 | SRH/Auger loop L2 |
|---:|---:|---:|---:|
| -19.7 V | 11268.75 | 2225.18 | 4.9961e-4 |
| -19.8 V | 11368.10 | 2225.97 | 5.0320e-4 |

At `-19.7 V`, transport is dominated by hole-from-hole (`8096.11`) and
electron-from-electron (`6729.26`) carrier blocks. Its largest entry is
`-3286.22`, from hole QFP node 7 `(0.75, 0.25)` to hole residual node 4
`(0.5, 0.25)`. At `-19.8 V`, the same support and sign remain, with value
`-3309.35`.

The avalanche loop is dominated by the electron-input paths. Its largest
entry is positive and remains on electron QFP node 13 `(1.25, 0.25)` to
electron residual node 10 `(1.0, 0.25)`: `805.53/802.47` at
`-19.7/-19.8 V`.

The transport loop has nearly balanced positive and negative L1 support
(`45688/45062` at `-19.7 V`), while avalanche is predominantly positive
(`11948` positive versus `2712` negative). Direction reversal therefore
cannot be inferred from a single matrix-entry sign; it is a global
carrier-coupled solve effect.

## Derivative verification

| Bias | `B` directional FD relative error | `C` directional FD relative error |
|---:|---:|---:|
| -19.7 V | 1.3073e-6 | 1.5859e-6 |
| -19.8 V | 9.0034e-7 | 2.2356e-6 |

Both checks pass the predeclared `1e-4` gate with a fixed `1e-7`
double-symmetric step. This rejects a gross sign, scale, or analytic
cross-derivative implementation error along the actual loop directions.
It does not prove every individual Jacobian entry or every nonsmooth branch.

## Scaling and conditioning

| Bias | raw `D` resolved condition | equilibrated `D` | equilibrated `S` |
|---:|---:|---:|---:|
| -19.7 V | 4.1247e13 | 2.9324e6 | 135.18 |
| -19.8 V | 4.1135e13 | 2.7976e6 | 137.60 |

`A` changes from a raw condition near `378` to an equilibrated condition near
`1.19`, showing a removable row/column scaling effect. The carrier-only `D`
block remains strongly ill-conditioned after equilibration. Closing the
Poisson-QFP loop reduces the equilibrated Schur condition to about `136`.

Consequently, the adverse QFP direction is not explained by the Schur system
being less numerically conditioned than the carrier-only system. The loop
regularizes the linear algebra while rotating the update away from the
Sentaurus target. The independent carrier direction remains evidence only:
its usable target-projection scale is about `2.7e-5/2.8e-5`.

## Closure and determinism

- `K` component closure is exact to exported precision at both biases.
- Schur/full-step agreement remains below `7.32e-13 V`.
- Schur relative closure remains below `4.60e-14`.
- Boundary QFP targets remain exactly zero at six contact nodes per bias.
- Independent roots have identical node, Jacobian-block, and Schur-loop CSV
  hashes at both biases.

Final Schur-loop hashes:

| Bias | Schur-loop CSV SHA-256 |
|---:|---|
| -19.7 V | `cdd875906efd403442a999ce8b8e3877c849f3289f1989da2f9da5732aaccbf6` |
| -19.8 V | `dabd7c4b2093df204aba78fa4acaafebff209936656a792eef6268d66d03af08` |

Generated scorecard:
`build-release/pn2d-task7-schur-source-final-scorecard-20260730/acceptance.json`.

## Decision

The typed cross-bias result is
`transport_and_avalanche_independently_sustain_reversal`.

The completed checks exclude continuation, step caps, direct contact rows,
SRH/Auger, gross `B/C` directional derivative errors, and a worse-conditioned
Schur system as primary explanations. They localize the remaining discrepancy
to how the production transport and avalanche cross-feedback paths select a
Newton direction relative to the proprietary target state.

Task 8 remains prohibited. Before proposing a correction, the next review
must compare the analytical ownership and sign conventions of the two
interior paths against the governing residual equations and an independent
reference definition. Removing a physical term, deleting a cross block, or
using the severely damped carrier-only direction is not authorized.
