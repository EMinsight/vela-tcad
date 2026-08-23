# SLOT-LDMOS strict IALMob BVDS continuation progress (2026-08-22)

## Outcome

The requested `IALMob off/on -> Id >= 1e-7 A/um -> interpolate inner BVDS`
run is not complete.  The direct-bordered physical load line is validated, but
the production `local_ad` device block still limits continuation to sub-nA
current increments near the folded avalanche branch.  No BVDS or IALMob shift
is reported until both threshold brackets have converged.

## Reproducible progress

| Item | Result |
|---|---|
| IALMob-off physical-R restart | Converged through inner voltage 15.829142186578611 V and drain current 1.4636193441414995e-9 A/um |
| IALMob-off direct-current continuation | Reached inner voltage 15.856737161516595 V and drain current 3.32361934414161e-9 A/um |
| Fold evidence | The accepted inner voltage moved 15.885851 -> 15.880551 -> 15.856737 V while current increased, confirming a negative-differential/folded branch |
| IALMob-on secant bootstrap | Fixed-inner bootstrap converged in 17 Newton iterations at 0.8078552725248964 V and 5.9253738367672866e-11 A/um |
| Direct-current scalar row | Added `coupled_voltage_coefficient`; `0` with `R=1` solves `I - I_target = 0`, while the default `1` preserves the physical series-resistor equation |
| Unit verification | `test_coupled_load_line.exe`: all 6 test cases and 41 assertions passed |
| Remaining threshold distance | off is still below the 1e-7 A/um criterion by about 30x; on production continuation has not started |

## Deterministic blocker

At an off target increment of 2e-11 A/um, the new direct-current row converged
to a scalar residual of 3.1e-23 A/um.  After 80 iterations the device block was
1.7745e-6, just above the strict 1e-6 production tolerance.  Relaxing the
locator device tolerance and using a frozen avalanche-source Jacobian allowed
additional accepted points, but multi-nA increments still spent many minutes
in device-block line search and retry.  Therefore the scalar boundary matrix is
not the active error; the remaining limitation is nonlinear globalization of
the full local-AD avalanche device block.

## Inexact device forcing implementation and smoke result

The augmented load-line Newton solver now has an opt-in dynamic device-block
tolerance.  It keeps the strict equation tolerance while the load row is far
from closure, then interpolates geometrically toward a configured maximum only
inside the load activation window.  Convergence, merit scaling, and residual-
filter scaling all use the same effective tolerance.  The default is disabled,
so existing production decks retain their previous strict behavior.

The bounded IALMob-off smoke run used strict/local-AD physics, an equation
tolerance of `1e-6`, a maximum forced tolerance of `5e-6`, a load tolerance of
`1e-12`, eight coupled iterations, and no step retries.  Results:

| Check | Result |
|---|---|
| Restart point, `Id=3.32361934414161e-9 A/um` | Accepted at iteration 0 with `coupled_inexact_device_forcing`; load residual `-8.27e-24` |
| Trial point, `Id=3.34361934414161e-9 A/um` | Correctly rejected after the bounded iteration budget |
| Final trial load residual | `8.567658926171108e-14`, below `1e-12` |
| Final trial device residual | `4.153328284762107`, far above `5e-6` |
| Interpretation | The scalar/load block is closed. Inexact forcing prevents unnecessary polishing of a nearly converged restart, but it cannot accept or conceal the unresolved local-AD device-block step. |

The smoke output is under
`outputs/ialmob_ablation/direct_bordered_20260822_v5/ialmob_off_direct_current_locator_inexact_smoke`.
The next solver task is therefore a device-block Jacobian/globalization audit
at the rejected trial, not a further relaxation of the convergence tolerance.

## Files

- `include/vela/simulation/CoupledLoadLine.h`
- `include/vela/simulation/DCSweep.h`
- `src/simulation/DCSweep.cpp`
- `tests/test_coupled_load_line.cpp`
- `scripts/prepare_slot_ldmos_direct_ialmob_bvds.py`
- `scripts/prepare_slot_ldmos_ialmob_high_r_locator.py`
- `scripts/prepare_slot_ldmos_frozen_state_eval.py`

The latest accepted off checkpoint is:

`outputs/ialmob_ablation/direct_bordered_20260822_v5/ialmob_off_direct_current_locator/states/state_bias_15p856737.csv`

## Required next implementation

1. Preserve the accepted current-controlled tangent across restarts using
   explicit target-current metadata rather than reconstructing targets from
   filenames.
2. Audit the rejected trial's local-AD device Jacobian and residual-filter
   search direction; do not increase the forcing cap above `5e-6` to hide a
   residual of order one.
3. Continue off and on to converged 8e-8 and 1.2e-7 A/um points.
4. Reclose those four states with the physical R=1e12 ohm.um, strict local-AD
   Jacobian, equation tolerance 1e-6, and filter envelope 2.
5. Interpolate each inner BVDS at 1e-7 A/um and report the IALMob shift.
