# TransportModels Id–Vg strict SRH and deep-off validation

## Acceptance rules

- A point is numerically resolved only when `|Id| >= 10 * |four-terminal KCL residual|`.
- Resolved deep-off points require log-current error <= 0.15 dex.
- The silicon SRH-generation/substrate-hole closure error must be <= 1%.
- Unresolved points are not assigned a passing relative-current comparison.

## Deep-off results

| Model | Vg (V) | Sentaurus Id (A/um) | Vela Id (A/um) | log error (dex) | SRH/substrate closure | Id/KCL | Status | Acceptance |
|---|---:|---:|---:|---:|---:|---:|---|---|
| DD | -1 | 1.635e-15 | 6.621e-24 | 8.393 | 0.0001959 | 1.212e-06 | numerically_unresolved | not_accepted |
| DD | -0.84 | 1.622e-15 | 6.621e-24 | 8.389 | 0.000207 | 5.549e-07 | numerically_unresolved | not_accepted |
| DD | -0.68 | 1.739e-15 | 6.62e-24 | 8.419 | 0.0003036 | 3.231e-08 | numerically_unresolved | not_accepted |
| DD | -0.52 | 1.616e-14 | 2.145e-14 | 0.1232 | 0.0003949 | 68.41 | resolved | pass |
| DG | -1 | 1.637e-15 | 4.471e-21 | 5.564 | 0.001008 | 0.0008116 | numerically_unresolved | not_accepted |
| DG | -0.84 | 1.638e-15 | 4.471e-21 | 5.564 | 0.001043 | 0.0003274 | numerically_unresolved | not_accepted |
| DG | -0.68 | 2.163e-15 | 3.274e-16 | 0.8199 | 0.001145 | 1.812 | numerically_unresolved | not_accepted |
| DG | -0.52 | 4.282e-14 | 3.95e-14 | 0.03508 | 0.001174 | 132.3 | resolved | pass |

## DD direction check

| Vg (V) | Forward Id (A/um) | Reverse Id (A/um) | Relative difference |
|---:|---:|---:|---:|
| 2.04 | 0.001662 | 0.001662 | 1.21e-10 |
| 1.88 | 0.001541 | 0.001541 | 1.544e-10 |
| 1.72 | 0.001411 | 0.001411 | 1.612e-10 |

The DD forward sweep completed to 2.2 V. The same-process reverse sweep matched the overlapping 2.04, 1.88 and 1.72 V points, then stopped near 1.623 V at the strict residual floor.

The DG sweep resolved the requested deep-off points through -0.52 V and was bridged to -0.40 V. It then stopped at -0.3987109375 V with residual 1.4374e-10 after adaptive step reduction; later DG biases are therefore not claimed as completed.
