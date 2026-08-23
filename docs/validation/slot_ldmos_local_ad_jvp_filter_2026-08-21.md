# SLOT-LDMOS local-AD, JVP, and coupled-filter validation (2026-08-21)

## Outcome

The three requested nonlinear-solver tasks are implemented:

1. `newton_jvp_probe` can select an exact zero- or one-based node id and an
   optional number of adjacent-cell rings. It can read a native Vela state CSV
   directly through `state_file`.
2. The canonical legacy `triangle_gss_gradqf_truncated` avalanche source has a
   nine-variable local forward-AD path selected by
   `impact_ionization.source_jacobian: local_ad`.
3. The coupled-load-line residual filter applies independent device and
   circuit envelopes. `residual_filter` no longer falls through to the global
   maximum-merit acceptance rule.

## Saved-state JVP audit

- State: `state_bias_15p720951.csv`
- Inner drain voltage: `15.720950570922257 V`
- Hotspot: zero-based node `10236`
- Perturbations: `psi`, `phin`, and `phip`
- Step sizes: `1e-4`, `1e-6`, and `1e-8 V`

| Jacobian mode | Maximum relative JVP error | Interpretation |
|---|---:|---|
| `local_ad` | `2.3974913934e-5` | Passes the `1e-3` target and `1e-4` preferred target |
| `frozen` | `3.2969961203e-2` | Missing avalanche derivatives are measurable at the hotspot |

At `1e-6 V`, the local-AD errors were `4.75e-8` for `psi`, `4.55e-7`
for `phin`, and `3.40e-8` for `phip`. The three step sizes form a stable
finite-difference consistency sweep rather than a single tuned comparison.

Generated evidence:

- `diagnostics/jvp/node10236_local_ad.csv`
- `diagnostics/jvp/node10236_frozen.csv`
- `simulation_jvp_node10236_local_ad.json`
- `simulation_jvp_node10236_frozen.json`

## Coupled-filter regression

The historical rejected transition was encoded as a unit regression:

- device residual: `8.2476318869e-7 -> 1.9599596828`
- load residual: `-1.6685485840 V -> -1.6620307948 V`
- trial damping: `1/256`

The new independent block filter rejects this transition. The prior shared
envelope used the very large normalized load residual to authorize the device
residual increase.

## Newly exposed continuation blocker

The saved state closes its load line at an outer voltage of approximately
`1187.82348632813 V`. A controlled `+1 V` outer step with `local_ad` and the
new filter does not yet converge: Eigen `SparseLU` fails while factoring the
fixed-inner-voltage device Jacobian, even after the outer controller reduces
the attempted step to `0.125 V`.

This failure is downstream of the three corrections. The present coupled
solver uses Schur elimination and therefore requires the device block `J` to
be invertible before it can form the circuit correction. A mathematically
correct full avalanche Jacobian can expose the fold singularity that the
frozen source Jacobian hid. The next development gate is therefore a direct
bordered solve (or equivalent robust pseudo-arclength formulation) for the
augmented device/circuit matrix, with explicit terminal-current derivatives.

## Tests

- `test_newton_solver.exe "[triangle_gss]"`: pass, 8 assertions.
- `test_coupled_load_line.exe`: pass, 30 assertions.
- `test_dc_sweep.exe "[coupled_newton]"`: pass, 42 assertions.
- Slot-LDMOS BVDS, external IALMob, and JVP preparation regression tests:
  pass, 5 tests.
