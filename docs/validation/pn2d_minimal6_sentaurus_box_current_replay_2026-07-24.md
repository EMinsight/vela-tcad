# PN2D Minimal6 Sentaurus box-current replay

Date: 2026-07-24

Status: `valid_reconstruction`

## Scope

This experiment reconstructs the Sentaurus O-2018.06-SP2 continuity-equation
edge current from documented box-method coefficients, endpoint state, and
element mobility. The result is an independently closed operator replay, not
a native edge-current dataset.

The initial state is Minimal6 `mirror` at -1 V. No Vela production formula
was modified.

## Documented discrete operator

For an oriented edge from vertex i to vertex j and adjacent element e,
Sentaurus uses:

```text
kappa_ij^e = d_ij^e / l_ij
B(x) = x / (exp(x)-1)
```

The Boltzmann endpoint state gives the following effective-potential
differences:

```text
x_n = ln(n_i/n_j) + (eQFP_i-eQFP_j)/Vt
x_p = ln(p_i/p_j) - (hQFP_i-hQFP_j)/Vt
```

The reconstructed element-edge currents are:

```text
F_n^e = q Vt kappa^e mu_n^e
        [n_i B(x_n) - n_j B(-x_n)]

F_p^e = q Vt kappa^e mu_p^e
        [p_j B(-x_p) - p_i B(x_p)]
```

The element contributions are summed for each global edge. The raw formula
has units A/cm in two dimensions and is multiplied by `1e-4` to obtain
A/um. The final replay uses the legacy Sentaurus elementary-charge constant
`q = 1.6021918e-19 C`. The modern exact SI value was retained as a control.

## Sentaurus self-closure

The reconstructed edge currents were summed with the exact runtime edge
orientation and compared with the high-precision Current file.

| Component | Reference, A/um | Reconstructed, A/um | Relative error |
|---|---:|---:|---:|
| Anode electron | -8.86921390144714e-18 | -8.86921291327539e-18 | 1.11416e-7 |
| Cathode electron | 1.29789524916234e-16 | 1.29789527054439e-16 | 1.64744e-8 |
| Anode hole | -1.25822279959092e-16 | -1.25822281777714e-16 | 1.44539e-8 |
| Cathode hole | 4.90196894430461e-18 | 4.90196841184066e-18 | 1.08622e-7 |

The maximum contact-component relative error is `1.11416e-7`.

The total electron-plus-hole KCL residuals at the two internal vertices are:

| Vertex | Residual, A/um |
|---:|---:|
| 1 | -3.49227e-25 |
| 5 | -4.26064e-25 |

The maximum KCL residual relative to terminal total current is
`3.16326e-9`.

Edges 1 and 6 have exactly zero accumulated box coefficient and exactly zero
reconstructed current for both carriers.

## Control branches

| Candidate | Contact relative-error range | Interpretation |
|---|---:|---|
| QFP plus for electrons, QFP minus for holes | 1.45e-8 to 1.11e-7 | Accepted discrete sign convention |
| Electrostatic-potential-only control | 0.4838% to 0.5188% | Close but not exact; omits endpoint state terms |
| Opposite QFP signs | 13x to 555x or worse | Rejected |
| Native element-vector projection | about 1.41e6 relative | Rejected as a box-edge current proxy |

The native element-vector projection differs from terminal-consistent box
current by approximately `6.15 dex`. This explains the earlier apparent
six-decade density/current inconsistency: element current vectors and box
continuity edge currents are different discrete observables.

## Initial same-edge Vela comparison

The reconstructed Sentaurus currents were mapped to the same nine canonical
edges as the authoritative Vela fixed-state `mirror/-1 V` replay. Both sides
use the Sentaurus box coefficient and one-micrometer depth convention.

Five edges per carrier are nonzero on both sides:

| Sentaurus edge | Vela edge | Carrier | Absolute error, dex | Sign agreement |
|---:|---:|---|---:|---:|
| 2 | 1 | electron | 1.42530 | yes |
| 2 | 1 | hole | 3.70136 | yes |
| 3 | 4 | electron | 1.39074 | yes |
| 3 | 4 | hole | 3.68905 | yes |
| 4 | 3 | electron | 5.11407 | no |
| 4 | 3 | hole | 5.20257 | no |
| 5 | 6 | electron | 3.68303 | yes |
| 5 | 6 | hole | 1.24088 | yes |
| 7 | 7 | electron | 3.68066 | yes |
| 7 | 7 | hole | 1.26703 | yes |

| Carrier | Valid edges | Median error, dex | Maximum error, dex | Sign agreement |
|---|---:|---:|---:|---:|
| electron | 5 | 3.68066 | 5.11407 | 80% |
| hole | 5 | 3.68905 | 5.20257 | 80% |

The majority-carrier horizontal transport retains the approximately
3.69-dex Vela/Sentaurus discrepancy. The central vertical edge has the only
sign disagreement and the largest magnitude gap. Minority-carrier contact
edges differ by approximately 1.24-1.43 dex.

This is the first comparison in this audit that uses identical directed-edge
topology and box geometry rather than a projected Sentaurus node or element
current.

## Decision

The Sentaurus box-current reconstruction is accepted as a terminal- and
KCL-closed diagnostic reference. It supersedes node-current projection and
element-vector projection for subsequent edge-current localization.

It does not yet authorize a Vela formula change. The remaining cross-solver
gap includes different solved QFP states, carrier densities, and mobility.
The next discriminating experiment is staged replacement on this exact
box-edge reference:

1. replace Sentaurus QFP only;
2. recompute density from potential, QFP, and effective intrinsic density;
3. replace element mobility;
4. replace box geometry;
5. compare the final Vela SG operator with the closed Sentaurus replay.

## Evidence

Deterministic evidence root:

`build-release/pn2d-minimal6-sentaurus-box-current-replay-20260724-a`

The root contains raw Sentaurus output, the high-precision Current file,
sealed Vela input, edge tables, closure summaries, scripts, hashes, and an
independent verifier. The independent verifier passed all contact, KCL,
zero-edge, and cross-solver-summary checks without importing the main
analyzer.
