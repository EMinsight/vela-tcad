# PN2D reverse mobility-basis qualification

Date: 2026-07-27

## Decision

`cell_reconstructed_total_impurity` does not pass the reverse-BV production
qualification yet. Do not change the global mobility default.

The candidate itself is numerically consistent and reduces the triangle GSS
source and source-Jacobian magnitude by about 1%. However, after the 2-D
charge-area correction, both the unchanged `net_doping` control and the
candidate enter the same nonphysical avalanche branch. The candidate therefore
cannot be judged against Sentaurus until the current production avalanche
baseline is repaired or requalified.

## Matched simulations

Four self-consistent simulations were run on the same coarse7x3 mesh and
production physics:

- `net_doping`, avalanche on
- `net_doping`, avalanche off
- `cell_reconstructed_total_impurity`, avalanche on
- `cell_reconstructed_total_impurity`, avalanche off

All simulations used 0.05 V internal continuation from 0 V to -20 V. Integer
bias points were extracted for comparison. This fine continuation was used
after a 1 V continuation control entered the same high-current branch.

## Reverse IV and breakdown gain

The gain criterion is `abs(I_avalanche_on) / abs(I_avalanche_off)`. Threshold
voltages are log-interpolated between integer bias points.

| Configuration | Gain 1.5 | Gain 2 | Gain 5 | Gain 20 |
|---|---:|---:|---:|---:|
| Sentaurus | 15.892 V | 17.367 V | 19.386 V | not reached |
| Vela `net_doping` | 3.042 V | 3.071 V | 3.165 V | 3.307 V |
| Vela cell-reconstructed | 3.042 V | 3.072 V | 3.166 V | 3.309 V |

The Vela avalanche-on branch rises to about `2.7e-9 A/um` near -15 V and
then falls to about `2.4e-14 A/um` at -20 V. This non-monotonic branch is not
comparable to the Sentaurus BV curve. At -20 V, the cell-reconstructed current
is still 55.4 times the Sentaurus current.

Changing the mobility basis shifts the gain-2 threshold by only 0.00048 V.
It does not correct the premature avalanche feedback.

## Avalanche source comparison

At the exact integer biases from -15 V to -20 V, the reconstructed basis
reduces the integrated triangle GSS source by approximately 0.97% to 1.29%
relative to `net_doping`. The terminal avalanche-on current changes by a
similar amount.

The weighted source centroid changes by less than `3.5e-5 um` in x and about
`1.5e-4 um` in y. The reported peak cell changes from cell 12/17 to cell 18;
the nearly unchanged centroid and source magnitude indicate an equivalent
local support selection rather than a material hotspot relocation.

## Source Jacobian gate

A new compensated-junction test evaluates the production triangle-GSS source
with both mobility bases. It subtracts an avalanche-off assembler, then
compares the assembled source Jacobian with an independent external central
finite-difference Jacobian.

| Basis | Assembled source-J norm | FD source-J norm | Relative difference |
|---|---:|---:|---:|
| `net_doping` | 2.481133999e16 | 2.481133936e16 | 4.0655e-6 |
| cell-reconstructed | 2.457526037e16 | 2.457525975e16 | 4.0610e-6 |

Both pass the `2e-4` gate. The candidate lowers the Jacobian norm by about
0.95%, consistent with its source reduction. There is no evidence of a
candidate-specific Jacobian defect.

## Interpretation and next action

The mobility-basis candidate is not the cause of the failed reverse curve:
the unchanged `net_doping` control fails in the same way. The current reverse
production baseline predates the 2-D charge-area correction, so its old curve
cannot be reused as a current-code control.

The next P0 task is to isolate the post-charge-area avalanche feedback by
comparing:

1. avalanche source integral versus continuity transport divergence;
2. source scaling in physical and unit-scaled coordinates;
3. fixed-state source evaluation on the post-fix reverse states;
4. triangle-GSS versus element-edge local-AD source paths.

Only after the `net_doping` production control again follows the Sentaurus
low-gain branch through roughly 17 V should the cell-reconstructed mobility
basis be requalified for the global default.

## Artifacts

- `build-release/pn2d-reverse-mobility-basis-qualification-fine-20260727/reverse_iv_gain_comparison.csv`
- `build-release/pn2d-reverse-mobility-basis-qualification-fine-20260727/breakdown_gain_thresholds.csv`
- `build-release/pn2d-reverse-mobility-basis-qualification-fine-20260727/avalanche_source_summary.csv`
- `build-release/pn2d-reverse-mobility-basis-qualification-fine-20260727/reverse_iv_gain_comparison.png`
- `tests/test_impact_ionization.cpp`
