# SingleDevice Eq. 231 Si/SiO2 affine reaction closure (2026-08-14)

## Outcome

The Si/SiO2 interface reaction discrepancy at nodes 848 and 2630 is closed
for both imported endpoint states. The Formula-0 fitted edge operator is not
the source of the error. The missing term is the constant part of the
algebraically eliminated material-side reaction trace, so the interface
reaction must be affine rather than a single control-volume multiplier.

## Direct Sentaurus probe

At linear-state node 848, eQuantumPotential was perturbed at the row and all
six one-ring nodes. Positive and negative `1e-5 V` perturbations were run in
Sentaurus O-2018.06 and the initial NewtonPlot residual was used for central
differences. Four edges have nonzero Formula-0 support; their Vela/Sentaurus
Jacobian ratios are:

- `-295.2907691`;
- `-295.2878666`;
- `-295.3029320`;
- `-295.2907666`.

A least-squares common normalization is `-295.3029099`. This proves that the
Vela edge flux has the same support and relative coefficients as Sentaurus.
The remaining row mismatch is diagonal reaction closure.

Combining the center-difference reaction Jacobian with all 69 free Si/SiO2
rows at each of the linear and saturation states gives one cross-state affine
trace:

| material side | Lambda slope | offset |
|---|---:|---:|
| Silicon | `0.3613278292533479` | `-0.00020247747279261268 V` |
| SiO2 at Silicon | `2.6839079693374917` | `-0.0015039829729206406 V` |

The direct Jacobian is reproduced to `1.3e-7` in the Vela-scaled row. The
cross-state fitted fixed-row error is below `8e-4`.

## Implementation

The four existing region-side reaction weights now have matching optional
offsets. A material-side residual uses

`coefficient * (weight * Lambda + offset)`

and the analytic diagonal Jacobian uses `coefficient * weight`; the constant
offset correctly has zero derivative. All offsets default to zero. The
SingleDevice experimental `sentaurus_box` profile sets only the Silicon and
insulator-at-Silicon offsets, leaving PolySilicon controls unchanged.

## Fixed-state and endpoint results

| endpoint | node 848 old/new residual | node 2630 old/new residual | raw outer change |
|---|---:|---:|---:|
| linear | `-0.445529 / -0.000617` | `-0.446803 / -0.000606` | `0.310 mV` |
| saturation | `-0.429258 / -0.000760` | `-0.445876 / -0.000614` | `0.516 mV` |

The linear one-step endpoint is below the `0.5 mV` acceptance threshold. The
saturation endpoint is `0.016 mV` above it, but its limiting update is now at
pure-Silicon node 1816; the Si/SiO2 interface rows are already below `8e-4`.
The largest fixed-state residual has moved to the PolySilicon/SiO2 family
(nodes 2075, 2072, and neighbors), about `0.332` linear and `0.302`
saturation. The next localization target is therefore not the Si/SiO2
control-volume reaction.

## Verification

- density-gradient quantum tests: 22 cases, 93 assertions passed;
- Newton/configuration tests: 83 cases, 1165 assertions passed;
- manufactured interface solve exercises nonzero affine offsets and verifies
  restart consistency;
- `git diff --check`: passed.
