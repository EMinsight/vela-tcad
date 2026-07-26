# PN2D exact high-bias Sentaurus oracle

Date: 2026-07-26

Plan task: Task 14 of
`docs/superpowers/plans/2026-07-26-pn2d-high-bias-process-variable-jacobian-localization.md`

## Outcome

Typed outcome: `exact_high_bias_oracle_available`

The exact coarse7x3 Sentaurus high-bias oracle is accepted. It does not
authorize Task 15 or any Vela nonlinear candidate because Task 12 stopped with
`incomplete_analytic_derivative`.

## Exact paired matrix

The lattice is:

`-10, -18, -19, -19.5, -19.8, -19.9, -19.95, -20 V`.

Nine otherwise controlled branches were run on two independent roots:

1. production implicit avalanche;
2. explicit `GradQuasiFermi`;
3. explicit electric-field avalanche drive;
4. explicit `GradQuasiFermi` plus `UseQuasiFermi` at contacts;
5. explicit electric-field drive plus `UseQuasiFermi` at contacts;
6. `AvalDensGradQF`;
7. high-field-saturation-disabled electric-field control;
8. high-field-saturation-disabled QFP-gradient control; and
9. avalanche disabled.

All `9 x 8 x 2 = 144` exact branch/state observations passed. No failed,
interpolated, or nearest-state row is present.

Raw roots:

- `build-release/pn2d-task14-high-bias-oracle-20260726-a`
- `build-release/pn2d-task14-high-bias-oracle-20260726-b`

For every branch, normalized runtime process records and exact CurrentPlot rows
are identical between A and B. The paired matrix uses the same mesh and
parameter hashes. Branch decks are generated from one template and differ
only by their declared avalanche drive, contact selector, HFS removal,
`AvalDensGradQF`, or avalanche removal.

## Exact implicit-default current turn

| Bias (V) | abs Anode total current (A) | source integral (A/um) |
|---:|---:|---:|
| -10 | 3.179640402e-17 | 7.306902619e-19 |
| -18 | 7.365358341e-17 | 4.244591356e-17 |
| -19 | 1.125860992e-16 | 8.130929396e-17 |
| -19.5 | 1.564824923e-16 | 1.248770659e-16 |
| -19.8 | 2.320220811e-16 | 1.958965435e-16 |
| -19.9 | 3.027800117e-16 | 2.608070677e-16 |
| -19.95 | 3.589097065e-16 | 3.126050609e-16 |
| -20 | 4.337466423e-16 | 3.823148999e-16 |

The turn is not a single discontinuous jump. Its log slope steepens
monotonically across the refined exact lattice:

| Interval (V) | current ratio | current log slope per V | source ratio | source log slope per V |
|---|---:|---:|---:|---:|
| -19 to -19.5 | 1.38989 | 0.65845 | 1.53583 | 0.85814 |
| -19.5 to -19.8 | 1.48274 | 1.31296 | 1.56872 | 1.50086 |
| -19.8 to -19.9 | 1.30496 | 2.66174 | 1.33135 | 2.86194 |
| -19.9 to -19.95 | 1.18538 | 3.40129 | 1.19861 | 3.62319 |
| -19.95 to -20 | 1.20851 | 3.78779 | 1.22300 | 4.02608 |

Thus the remembered coarse7x3 knee is present: the same voltage-normalized
current growth becomes much steeper as the state approaches -20 V.

## First accelerating process stage

For the implicit branch, exact -19 to -20 V growth is:

| Dependency quantity | -19 V | -20 V | ratio |
|---|---:|---:|---:|
| maximum electric field (V/cm) | 3.90538e5 | 4.01378e5 | 1.02776 |
| maximum electron QFP gradient (V/cm) | 3.88856e5 | 4.00165e5 | 1.02908 |
| maximum electron mobility (cm2/V/s) | 727.054 | 727.054 | 1.00000 |
| maximum electron velocity (cm/s) | 2.14389e6 | 4.54900e6 | 2.12184 |
| maximum element electron current density (A/cm2) | 3.20171e-2 | 1.23604e-1 | 3.86056 |
| maximum electron alpha (cm-1) | 2.68794e4 | 3.07343e4 | 1.14342 |
| maximum total generation (cm-3 s-1) | 2.02759e15 | 9.53368e15 | 4.70196 |
| source integral | 8.13093e-17 | 3.82315e-16 | 4.70198 |
| Anode total current | 1.12586e-16 | 4.33747e-16 | 3.85258 |

The first strong acceleration in the exported dependency order is therefore
the velocity/current layer. Maximum field and QFP-gradient change by about
3%, the maximum mobility does not change, and maximum alpha changes by about
14%; velocity doubles and element current increases by about 3.86x before
generation/source increases by about 4.70x and terminal current by about
3.85x.

This is a stage-order observation, not a proof of a proprietary local
constitutive formula. Task 15 would still be required to distinguish density,
velocity, mobility-support, current-support, and source-mapping effects at
matching local support.

Native electron, hole, and mean ionization-integral fields are present but
are exactly zero on this accepted run. They remain native-zero observations
and are not converted to finite log errors.

## Physics controls

At -20 V:

| Branch | source integral (A/um) | abs Anode total current (A) |
|---|---:|---:|
| implicit default | 3.82315e-16 | 4.33747e-16 |
| explicit GradQF | 3.82315e-16 | 4.33747e-16 |
| explicit electric field | 2.92823e-15 | 2.97966e-15 |
| GradQF + contact QFP selector | 3.82166e-16 | 4.33587e-16 |
| electric field + contact QFP selector | 2.92722e-15 | 2.97864e-15 |
| AvalDensGradQF | 9.30155e-17 | 6.06956e-17 |
| HFS disabled, electric-field drive | 2.53367e-15 | 2.58517e-15 |
| HFS disabled, QFP-gradient drive | 3.15916e-16 | 3.67419e-16 |
| avalanche disabled | 0 | 5.14294e-17 |

Implicit default and explicit `GradQuasiFermi` are exactly identical across
the complete summarized lattice. The global QFP-gradient default is therefore
retained. Electric-field drive is a materially different, much stronger
branch and is not substituted. The contact selector changes the -20 V
results by less than 0.04% and is not the knee cause. Removing avalanche
reduces the -20 V terminal current to about 11.9% of the implicit result,
which directly confirms that avalanche feedback is essential to the turn.

## Preserved observation boundary

- Current vectors are native node/element vectors.
- `ReadCoefficient` and `ReadMeasure` are native documented mesh operators.
- A directed-edge carrier current reconstructed from state and coefficients
  remains `operator_replay`.
- Native directed-edge current and native element alpha remain unsupported.
- No ordinary Plot/CurrentPlot residual or Jacobian is claimed.
- No Van Overstraeten, mobility, field, alpha, geometry, current, source, or
  voltage scale was fitted.

## Generated analysis

The independent analyzer writes:

- `process_summary.csv`;
- `adjacent_growth.csv`; and
- `reproducibility.csv`

under
`build-release/pn2d-task14-high-bias-oracle-analysis-20260726`.

All nine reproducibility rows pass both normalized runtime-record equality and
CurrentPlot exact-row equality.
