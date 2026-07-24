# PN2D Minimal6 electric-field mobility candidate

Date: 2026-07-24

Status: complete with typed outcome `config_improves_but_misses_target`.

## Decision

The config-only candidate changes only:

`solver.mobility.high_field_driving_force:
quasi_fermi_gradient -> electric_field`.

The impact-ionization configuration, mobility formula, saturation velocity,
Scharfetter-Gummel formula, source model, continuation, and tolerances remain
unchanged.

Both topologies converge at every exact `-1..-20 V` checkpoint. The candidate
improves electron and hole QFP median and P95 errors on both mirror and sketch,
but it misses every frozen QFP, density, current, terminal-current, and impact
target. It is therefore retained as a comparison configuration and is not
adopted as the production default.

## Deterministic sweep evidence

Candidate sweep roots:

- `build-release/pn2d-minimal6-electric-field-sweep-20260724-a`
- `build-release/pn2d-minimal6-electric-field-sweep-20260724-b`

Each root contains 40 accepted Vela checkpoints, 42 accepted Sentaurus
checkpoints including the two 0 V starts, and zero rejected transitions.

Final production-source Phase F roots:

- baseline:
  `build-release/pn2d-minimal6-phase-f-triangle-source-20260724-a/b`;
- electric-field candidate:
  `build-release/pn2d-minimal6-electric-field-phase-f-triangle-source-20260724-a/b`.

The independent verifier passes both A/B pairs with 40 states, 80 internal-node
rows, 320 carrier-element mobility rows, 400 directed carrier-edge rows, and
40 terminal/source rows.

## Dependency-chain comparison

| Quantity | QFP baseline median | Electric candidate median | Candidate result |
|---|---:|---:|---|
| electrostatic potential, V | `5.36593e-12` | `5.36593e-12` | unchanged, pass |
| electron QFP, V | `0.0548737` | `0.0514769` | improved, target missed |
| hole QFP, V | `0.0620671` | `0.0591186` | improved, target missed |
| electron density, dex | `0.921833` | `0.864770` | improved, target missed |
| hole density, dex | `1.042676` | `0.993144` | improved, target missed |
| electron element mobility, dex | `0.465174` | `0.465257` | essentially unchanged |
| hole element mobility, dex | `0.292739` | `0.292733` | essentially unchanged |
| electron directed current, dex | `0.823474` | `0.823152` | small improvement |
| hole directed current, dex | `0.939922` | `0.939636` | small improvement |
| terminal current, dex | `0.647139` | `0.646672` | small improvement |
| production triangle impact source, dex | `10.757025` | `10.764920` | slightly worse |

All 80 internal electron-QFP samples and all 80 internal hole-QFP samples
improve in paired comparisons. Mirror and sketch give the same conclusion
after the verified cell permutation.

Directed-current sign agreement remains 80 percent because the 80 central
`1-5` carrier-edge samples retain the previously classified near-zero
reference sign mismatch. The candidate does not create a new sign class.

## Source metric correction

Two diagnostic source-contract defects were found while executing this task:

1. Sentaurus source integration treated package mesh coordinates stored in
   micrometers as meters, inflating area and source by `1e12`.
2. Phase F used the sweep global-edge source proxy instead of the production
   `triangle_gss_gradqf_truncated` source.

The final roots above honor `mesh.coordinate_unit` and integrate the raw
triangle operator-audit source. Earlier Phase F roots reporting approximately
`12.9 dex` or `0.895 dex` impact errors are excluded from the final impact
decision.

No production C++ formula is changed by this candidate experiment.
