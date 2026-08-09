# BV Jacobian fixed-pattern assembly result

## Outcome

The three-step Jacobian work is complete:

1. Phased instrumentation established that edge-physics derivatives consume about 91% of Jacobian time, node sources about 9%, and triplet finalization less than 0.25%.
2. Mesh adjacency and a boundary-aware compressed sparse pattern are cached, so Newton state changes no longer invalidate the symbolic pattern.
3. Edge, node, and cell scatter offsets are precomputed and contributions are accumulated directly through `SparseMatrix::valuePtr()` after zero-valued pattern initialization.

Dirichlet treatment is unchanged. Constrained equations remain explicit unit-diagonal rows; unknowns are not eliminated and no large diagonal penalty is introduced.

## Structural result

| Scenario | Jacobians | Pattern builds | Pattern changes | Analyze calls | Analyze cache hits | Structural nnz |
|---|---:|---:|---:|---:|---:|---:|
| 6.08709 -> 6.09959 V | 192 | 1 | 0 | 5 | 191 | 103,229 |
| Voltage-to-Current final | 280 | 2 | 0 | 2 | 278 | 103,229 |
| External resistor 1206 V | 663 | 4 | 0 | 5 | 659 | 103,229 |

The cell/material-aware stencil tightening reduced structural nnz from the step-2 superset of 133,943 to 103,229. Against the numerical nnz of about 82--83k, structural overhead is now about 24--26%, rather than about 63%.

Triplet storage is no longer used during Newton assembly (`triplet_capacity = 0`). The former finalization pass fell from 1.14--3.65 seconds to 8--33 microseconds over an entire benchmark, a reduction greater than 99.999% in all three scenarios. SparseLU symbolic analyze calls fell by 97.3--99.3% relative to the state-dependent baseline.

## Numerical gates

| Gate | Candidate |
|---|---:|
| Voltage-to-Current BV | 6.395904174606575 V |
| Voltage-to-Current boundary residual | 6.169747e-12 A/um |
| External-resistor BV | 6.395887865540998 V |
| External-resistor load-line residual | 9.136147e-8 V |
| V2I electron/hole continuity | 1.216560e-9 / 1.308305e-9 |
| External electron/hole continuity | 3.115851e-10 / 3.237868e-10 |

The V2I field comparison also remains at the prior gate: potential correlation 0.999996, populated-carrier electron and hole quasi-Fermi correlations 0.999996 and 0.999990, and high-field edge correlation 0.992835.

No mobility, avalanche, BTBT, empirical-current-scaling, convergence-tolerance, or boundary-equation parameter was changed.

## Performance interpretation

Single-run wall/Jacobian timing moved in inconsistent directions: high-field transition improved about 12%, V2I regressed about 23%, and external resistor regressed about 2%. These runs therefore do not demonstrate a stable end-to-end speedup on the Intel N150 host. The deterministic gains are narrower: fixed symbolic structure, 97--99% fewer analyze calls, zero triplet allocation, and elimination of triplet merging.

The next justified optimization target is edge-physics derivative evaluation, which still accounts for about 91% of Jacobian assembly. It should be optimized only with numerical-equivalence tests because it contains the transport and avalanche derivative paths.

Machine-readable evidence is in `docs/validation/bv_jacobian_assembly_step3_2026-08-09.json`; raw profiles are under `build-bv-performance/jacobian-step3-fixed-scatter-tight`.
