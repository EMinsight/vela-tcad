# TransportModels DG self-consistent discretization validation

Endpoint: `Vg=1 V`, `Vd=2 V`; corrected material contract and neutral interface.
Sentaurus reference drain current: `0.000705525753105 A/um`.

| Operator | Converged | Id (A/um) | Absolute relative error | Iterations |
|---|---:|---:|---:|---:|
| P1 direct control | False | n/a | n/a | n/a |
| Sentaurus box contender | True | 0.000712295708665 | 0.9596% | 26 |

## Decision

- Lowest converged endpoint-current error: **Sentaurus box contender**.
- A fixed-state residual improvement is not accepted as production evidence unless
  this self-consistent check also converges and preserves terminal-current accuracy.
