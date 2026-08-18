# BVmethods non-transient freeze

Status: **PASS_WITH_CONTINUATION_PENDING**

Transient node 8 is explicitly outside this scope.

| Method | Node | Rows | Sentaurus BV | Reference frozen | Vela mapping |
| --- | ---: | ---: | ---: | --- | --- |
| ABA_poisson | 3 | 24 | 5.30552563 V | yes | operator reference |
| ABA_coupled | 4 | 28 | 6.37749428 V | yes | path/current IIC |
| resistor | 5 | 63 | 6.37979164 V | yes | external resistor |
| voltage2current | 6 | 89 | 6.3831842 V | yes | voltage to current |
| continuation | 7 | 59 | 6.38372717 V | yes | pseudo-arclength |

## Vela acceptance

- Path IIC: `pass`.
- Current IIC: `pass_within_3_percent`.
- External resistor: `PASS`.
- Voltage to current: `PASS`.
- Continuation: `pending`.

The Sentaurus continuation reference is frozen, but Vela continuation is not claimed complete until a converged NMOS arclength curve is supplied.

The bounded numerical trials are recorded in `continuation_diagnostic_20260818.json`; no physical parameter was changed.
