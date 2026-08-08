# BVmethods NMOS boundary-method implementation

## Frozen baseline

The implementation starts from local `main` commit `38ede5c`. The path-IIC and
current-IIC baseline in
`docs/validation/bvmethods_nmos_iic_final_validation_2026-08-06.md` was
regenerated before circuit work. Its seven focused Python unit tests passed.
No mobility coefficient, impact-ionization coefficient, or empirical current
scale was changed.

## Sentaurus external-resistor semantics

The official external-resistor deck is `full_raw/pp5_des.cmd`. It declares
`Resistor=1e7` on the drain. Direct checks of `n5_des.plt` close

```text
drain OuterVoltage = drain InnerVoltage
                   + 1e7 ohm*um * drain TotalCurrent [A/um]
```

to about `1e-12 V`. Drain current is positive on the reference branch.
`InnerVoltage` is the device voltage used for BV extraction; `OuterVoltage`
includes the resistor drop. At `|Id|=1e-4 A/um`, Sentaurus gives
`6.379791636301563 V`.

The voltage-to-current deck is `full_raw/pp6_des.cmd`. It ramps the drain to
`6 V`, executes `Set ("drain" mode current)`, and then ramps drain current.
At the same current threshold Sentaurus gives `6.38318420057198 V`.

## Implementation and tests

`BoundaryControl` provides the bracketed scalar closure used by both modes.
`DCSweep` adds mutually exclusive `external_circuit` and
`voltage_to_current` configurations and writes explicit inner/outer voltage,
target-current, resistance, residual, and evaluation-count columns. Boundary
evaluations are persisted as state CSV files. A restarted process restores the
best point and both sides of an existing bracket, then resumes from the secant
prediction instead of repeating the completed DD solves.

Both boundary modes reuse the configured continuation predictor for `psi`,
`phin`, and `phip`. The predictor uses the last two meaningful converged states
across voltage points and across current targets. If a predicted DD solve
fails, the same voltage is retried from the last constant warm state. Every
predicted state is still fully Newton-corrected and must pass the global
electron/hole continuity gate before it is accepted.

Focused Catch2 results:

- external-resistor equation, prediction, persistence, resume, and DCSweep
  integration: 77 assertions in 6 test cases;
- current-boundary root, state prediction, and voltage-to-current integration:
  17 assertions in 2 test cases;
- impact-ionization regression suite: 734 assertions in 52 test cases;
- all focused assertions passed on 2026-08-08.

The reproducible BVmethods driver is
`scripts/run_bvmethods_nmos_boundary_methods.py`.

## Self-consistent activation finding

A direct high-bias activation probe at `5.93977763593775 V` reduced the
combined Newton residual from `56.2994` to `6.0069e-10`, but stopped with
`carrier_row_convergence_line_search_rejected`: 120 carrier rows remained,
with maximum row ratio `0.1807705`. This was reproduced without the external
resistor loop, so it is an avalanche-feedback branch-initialization issue, not
a load-line defect.

The same final IIC operator converged at `0 V` from the accepted low-bias
state. Therefore production comparison builds the self-consistent avalanche
branch from low bias before entering the resistor or current-driven segment,
matching the Sentaurus solve order. The carrier-row ratio is retained as a
`report`-mode diagnostic for this branch. Enforced global electron/hole
continuity closure is the hard physical acceptance gate. This does not change
the self-consistent avalanche equations, mobility, or avalanche parameters.

That branch strategy accepted `0`, `0.01`, `0.05`, `0.1`, `0.2`, `0.35`,
`0.5`, `0.75`, and `1.0 V`. The runner discovers these per-point state files
and resumes from the highest completed nominal prebias point, so the remaining
path can continue from `1.0 V` without repeating the accepted points.
An attempted direct `1.0 -> 1.5 V` step remained in Newton iteration after
about `263 s` of CPU. The remaining prebias lattice therefore uses `0.1 V`
nominal spacing; this is a numerical continuation change only.
The low-bias handoff sweep also enables the existing secant state predictor for
`psi`, `phin`, and `phip`; every predicted point is still fully Newton-corrected
before it can become an accepted checkpoint.

The completed `0 -> 5.9 V` prebias contained 50 converged points. The first
external-resistor evaluation then rejected the already converged `5.9 V`
handoff solely through `carrier_row_convergence_line_search_rejected` when the
auxiliary row gate was changed back to `enforce`. The production comparison
therefore retains full global Newton convergence and full-strength avalanche
equations while recording the carrier-row gate in `report` mode. This is a
solver acceptance-policy change, not a mobility or avalanche-parameter fit.

The first load-line attempt completed ten device evaluations without yet
bracketing its first outer-voltage point because the scalar bracketing step was
limited to `0.025 V`; each high-field DD evaluation took roughly six minutes.
The runner now uses `0.1 V`, matching the accepted high-field prebias lattice,
while leaving the load-line residual tolerance unchanged.

The boundary workflow now treats the local carrier-row ratio as an audit and
uses enforced global electron/hole continuity closure as the hard conservation
gate. This separates near-zero row scaling sensitivity from physical current
conservation without relaxing the coupled Newton residual or changing any
transport or avalanche parameter.

## Final validation

The external-resistor run uses `1e7 ohm*um` and closes the Sentaurus load-line
sign convention. At `1e-4 A/um` it gives:

| Method | Vela BV (V) | Sentaurus BV (V) | Delta (V) | Relative error | 3% gate |
| --- | ---: | ---: | ---: | ---: | --- |
| External resistor | 6.395887865509238 | 6.379791636301563 | 0.016096229207674 | 0.252300% | PASS |
| Voltage to current | 6.395904174606575 | 6.383184200571980 | 0.012719974034595 | 0.199273% | PASS |

For the voltage-to-current branch, the final independently solved current
boundary point is `6.395904200106065 V`. Its measured current is
`1.0000000616974712e-4 A/um`, so the current-boundary residual is
`6.169747113378421e-12 A/um`. The global electron and hole continuity ratios
are `1.2165594664010111e-9` and `1.3083069725092697e-9`, respectively. The
reported BV is linearly extracted at exactly `1e-4 A/um`, yielding
`6.395904174606575 V`.

Machine-readable summaries are written under:

- `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/boundary_external_resistor_20260806/summary.json`
- `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/boundary_voltage_to_current_20260806/summary.json`

Both summaries record `pass_3_percent: true` and
`mobility_or_avalanche_parameter_fit: false`.

## Same-current Sentaurus field validation

An independent Sentaurus O-2018.06-SP2 current-boundary run was completed on
2026-08-08. It first ramped the drain voltage to `6.0 V`, changed the drain
electrode from voltage mode to current mode, and solved the fully coupled
Poisson/electron/hole system to exactly `1e-4 A/um`. The final imported contact
state is `6.384111661907364 V` and `1e-4 A/um`. This voltage is `0.014530%`
above the archived `6.383184200571980 V` Sentaurus result. The independently
solved Vela boundary point is `0.184717%` above this fresh Sentaurus result.

The terminal TDR contains 1909 semiconductor nodes that map exactly to the
Vela mesh by node ID. At these matching nodes, the potential median, 95th
percentile, and maximum absolute errors are `0.000764 V`, `0.015811 V`, and
`0.052573 V`; the spatial correlation is `0.9999957`. The median potential
offset is only `5.46e-6 V`, so a global gauge shift does not explain the
remaining high-field shoulder.

For carrier-populated nodes, defined independently for each carrier as a
Sentaurus density at least `1e-6` of its peak density:

| Quantity | Populated nodes | Median error | 95th-percentile error | Correlation |
| --- | ---: | ---: | ---: | ---: |
| Electron density | 1310 | 0.003353 dex | 0.078628 dex | 0.999810 |
| Hole density | 1112 | 0.023540 dex | 0.325500 dex | 0.996756 |
| Electron quasi-Fermi potential | 1310 | 0.002900 V | 0.016477 V | 0.999996 |
| Hole quasi-Fermi potential | 1112 | 0.005291 V | 0.026275 V | 0.999990 |

The absolute log-density errors correspond to median multiplicative mismatch
factors of `1.00775` for electrons and `1.05570` for holes. The
95th-percentile mismatch factors are `1.19847` and `2.11593`, respectively.
The larger hole tail is localized
relative to the overall distribution; its log-density correlation remains
`0.9968`.

Electric field is compared as the absolute projection along 5592 matching
semiconductor mesh edges. Sentaurus is averaged from its two nodal field
vectors and projected onto the edge; Vela is computed directly as
`abs(delta(psi))/edge_length`. Sentaurus and Vela peak projected fields are
`2.297503e8 V/m` and `2.574820e8 V/m`, so Vela is `12.0704%` higher at the
peak. On edges above 10% of the Sentaurus peak, the median and 95th-percentile
relative errors are `2.7639%` and `17.2860%`, with correlation `0.992835`.
Thus the high-field shape is strongly aligned, while a measurable peak-field
and hole-density shoulder remains. This is consistent with the previously
declared cross-simulator high-field state residual and was obtained without
mobility or avalanche-parameter fitting.

The comparison is regenerated by
`scripts/compare_bvmethods_nmos_boundary_fields.py`. Its machine-readable
output is
`build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_boundary_state_20260808/analysis/vela_sentaurus_field_comparison.json`.
