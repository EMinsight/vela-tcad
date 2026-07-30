# PN2D Task 7 Low-Current Root-Cause Audit

Date: 2026-07-30
Scope: observation only; no production default or collision-ionization physics
setting was changed.

## Typed outcome

`low_current_state_precision_floor_not_avalanche_operator_or_terminal_extractor`

The three low-current reverse intervals that blocked the predeclared Task 7
no-relocation gate are not evidence of a nonmonotone impact-ionization
coefficient or source. They occur at a terminal-current level of approximately
`3e-17 A/um`, where the solution is nearly in transport equilibrium and the
contact quasi-Fermi drop is approximately `3e-14 V`.

Two of the three intervals are already present in both avalanche-off and IIC.
The remaining avalanche-on-only interval moves when only Newton tolerances are
changed. The strict solve terminates with `stall_residual_floor`, which directly
identifies the numerical layer responsible for the interval pattern.

Task 8 and production-default changes remain unauthorized because this audit
does not retroactively change the original no-relocation acceptance gate.

## Experiment contract

The audit used the already sealed SG/Laux candidate:

- `current_approximation = element_edge_sg_gss_laux`;
- `source_mapping_mode = element_vertex_box_measure`;
- `driving_force = quasi_fermi_gradient`;
- `minimum_field_V_m = 0`;
- electron and hole driving-force reference densities equal to zero;
- the same 29-point Sentaurus exact bias lattice;
- the same `standard_0p05` continuation schedule.

Three independent controls were added:

1. duplicate self-consistent avalanche-on run;
2. SG-flux versus continuity-residual terminal-current extraction;
3. strict Newton tolerance run with `reltol=1e-12` and `abstol=1e-13`.

The controls are opt-in diagnostics. They do not alter the production default.

## Low-current evidence

| Bias (V) | Off current (A/um) | IIC current (A/um) | On current (A/um) | Max alpha (1/m) | Integrated qG proxy |
|---:|---:|---:|---:|---:|---:|
| -3 | -3.092773e-17 | -3.092773e-17 | -3.092773e-17 | 2.865282e-2 | 3.876163e-13 |
| -4 | -3.047925e-17 | -3.047925e-17 | -3.070629e-17 | 1.115797 | 1.741160e-11 |
| -5 | -3.048507e-17 | -3.048507e-17 | -3.048507e-17 | 11.916313 | 3.888583e-10 |
| -6 | -3.139930e-17 | -3.139930e-17 | -3.139930e-17 | 62.415966 | 4.709965e-9 |
| -7 | -3.049748e-17 | -3.049748e-17 | -3.049749e-17 | 211.934114 | 2.922258e-8 |

The collision-ionization coefficient and raw source both increase strictly with
reverse bias. There is no source kink at any reverse interval.

The hole drift and diffusion diagnostic components are individually
approximately `4e-11` to `1e-10 A/um`, while the net hole current is
approximately `3e-17 A/um`. Their equilibrium-balance ratio grows from
`2.68e6` to `6.40e6`. The production terminal current itself is computed from
the cancellation-free quasi-Fermi SG flux; the component ratio therefore
describes physical conditioning rather than a direct subtraction defect in the
terminal-current extractor.

## Branch and determinism evidence

- Avalanche-off and IIC state files are byte-identical at every audited bias.
- Their terminal currents are identical.
- Duplicate avalanche-on state files are byte-identical.
- No branch has a continuation retry at the audited exact targets.
- Standard-tolerance reverse intervals:
  - off/IIC: `-3 -> -4 V`, `-6 -> -7 V`;
  - on: `-3 -> -4 V`, `-4 -> -5 V`, `-6 -> -7 V`.

Therefore two intervals pre-exist without avalanche feedback, and the result is
deterministic rather than continuation-path random noise.

## Terminal-current extraction control

For both off and on branches over `-3` through `-7 V`:

- maximum relative difference between SG-flux and residual-derived terminal
  current: `4.650654e-15`;
- maximum absolute difference: `1.417484e-31 A/um`;
- maximum contact hole quasi-Fermi drop: `3.019807e-14 V`;
- diagnostic QF-floor current minus normal SG current: exactly zero.

The terminal-current extraction implementation is therefore excluded as the
cause of the reverse intervals at the measured precision.

## Strict Newton control

The strict run kept the physical operator and continuation schedule fixed and
changed only the Newton tolerances.

| Branch | Standard reverse intervals | Strict reverse intervals |
|---|---|---|
| off | `-3 -> -4`, `-6 -> -7` | `-3 -> -4`, `-5 -> -6` |
| on | `-3 -> -4`, `-4 -> -5`, `-6 -> -7` | `-3 -> -4` |

At the five exact targets, strict runs took 8-30 Newton iterations instead of
two. Every one terminated with `stall_residual_floor`; final residual norms
were `1.99e-12` to `1.05e-11`. The contact quasi-Fermi drop remained on the
same approximately `3e-14 V` plateau.

The movement and disappearance of intervals under a solver-tolerance-only
control proves that the low-current pattern is a state precision/stagnation
artifact. It is not a stable property of the impact-ionization operator.

## Sentaurus IIC semantics verified over SSH

The local Sentaurus training page
`/home/tcad/Desktop/Sentaurus_Training/sd/sd_11.html` states that IIC:

- solves the complete Poisson, electron, and hole transport system;
- computes ionization integrals and avalanche generation in postprocessing;
- uses `AvalPostProcessing` to exclude avalanche generation from the
  self-consistent carrier continuity equations;
- can record `AvalancheGeneration(Integrate(Semiconductor))`,
  `ImpactIonization(Maximum(Semiconductor Coordinates))`, and
  `ElectricField(Maximum(Semiconductor Coordinates))`.

The observed byte-identical off/IIC states are the expected implementation of
that contract.

## Open-source implementation comparison

The local repositories were audited at these immutable revisions:

| Code | Revision | Relevant behavior |
|---|---|---|
| Genius-TCAD-Open | `543da8452d5dfd33e6f8c457f962f6f670f0fce7` | Reconstructs cell current directions, uses directed SG edge currents, evaluates `G = alpha * abs(J) / q`, and partitions generated charge directionally on triangular elements. |
| Charon | `7cc38745625a6011ae3584ed111ec7ee74fb890e` | Uses vector current magnitudes, selectable QFP/current-aligned driving fields, optional density damping, and a minimum-field switch. |
| PISCES | `7f752a4ee279dad585d37e43a86287f185e89cdc` | Uses absolute directed side current, side area, field-dependent alpha, and analytic source derivatives before distributing the result to triangle nodes. |
| DEVSIM | `58a9a87083db00c6cadc0b4011c801db2cec5844` | Provides a user-defined finite-volume PDE framework; no packaged avalanche implementation was found in the local tree, so it is not a direct formula oracle. |

The common structure in Genius, Charon, PISCES, the Laux-Grossman formulation,
and Vela's SG/Laux candidate is:

`G = alpha_n(F_n) * |J_n| / q + alpha_p(F_p) * |J_p| / q`,

with multidimensional current support and a geometry-consistent distribution
of generated charge. This supports the earlier conclusion that Vela's former
midpoint-only current proxy had a support/semantic mismatch, while the complete
SG/Laux vector candidate is the physically justified correction.

The open-source codes also show optional low-field robustness controls:

- Genius disables and then smoothly ramps ionization below fractions of a
  critical field;
- Charon exposes a minimum field and density-damped driving force.

These are not evidence that the Sentaurus reference used either control.
Adding either after seeing the Task 7 result would be a new physics/numerics
axis and is prohibited by the present gate.

Primary references:

- Laux and Grossman,
  https://research.ibm.com/publications/a-general-control-volume-formulation-for-modeling-impact-ionization-in-semiconductor-transport--1
- Mauri et al., https://arxiv.org/abs/1412.3691
- Sandia Charon User Manual,
  https://www.sandia.gov/app/uploads/sites/106/2021/09/Charon_UserManual.pdf
- DEVSIM repository, https://github.com/devsim/devsim
- PISCES repository, https://github.com/ComputerWhisperer/pisces

## Root-cause statement

Two distinct issues must not be conflated:

1. **BV-active impact-ionization mismatch:** the dominant difference was the
   multidimensional current-density support and element-to-node source
   mapping. The opt-in complete element-edge SG/GSS-Laux candidate fixed the
   knee-region comparison from `11.400736` to `0.0146388 dex` RMSE and reduced
   breakdown-voltage error from `0.232 V` to `0.021 V`.
2. **Residual low-current nonmonotonicity:** this occurs near `3e-17 A/um`,
   is shared with avalanche-off/IIC, moves under Newton-tolerance-only
   controls, and ends in `stall_residual_floor`. It is a low-current state
   precision limitation, not a collision-ionization intermediate-quantity
   error.

Remaining collision-ionization comparisons should therefore stay focused on
the BV-active region and on native element current vectors, driving fields,
alpha, local generation, integral source, and source-to-node mapping. A
low-field current monotonicity defect must be handled as a separate nonlinear
precision task.

## Code recommendations

1. Keep the new terminal-current method comparison and Newton tolerance
   overrides diagnostic-only.
2. Preserve both raw acceptance results:
   `tradeoff_without_parity` for the original predeclared gate and
   `low_current_state_precision_floor_not_avalanche_operator_or_terminal_extractor`
   for the causal audit.
3. Add a prospective, reviewed acceptance contract that reports low-current
   precision-floor intervals separately from the BV-active curve/knee metrics;
   do not edit the old score retroactively.
4. For a separate solver task, add a low-current convergence observable based
   on contact QFP drop/current stationarity and carrier-block residuals.
5. Do not introduce a minimum field, density damping, source scale, or
   production-default change without a separately predeclared experiment.

## Artifacts

- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/acceptance.json`
- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/bias_diagnostics.csv`
- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/state_differences.csv`
- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/continuation_diagnostics.csv`
- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/terminal_current_method_compare.csv`
- `build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/tolerance_sensitivity.csv`
- `build-release/pn2d-task7-terminal-current-method-compare-20260730/`
- `build-release/pn2d-task7-strict-newton-low-current-control-20260730/`
