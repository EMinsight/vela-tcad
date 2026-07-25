# PN2D Minimal6 element-edge self-consistent candidate

Date: 2026-07-25

Status: complete with typed outcome
`config_improves_but_misses_target`.

## Candidate

The candidate changes the existing electric-field mobility comparison branch
by explicitly selecting the opt-in avalanche operator:

- `current_approximation: element_edge_sg_gss_laux`
- `driving_force: electric_field`
- `source_mapping_mode: element_vertex_box_measure`
- `quasi_fermi_gradient_discretization: edge_difference`

The mobility branch remains `masetti_field` with
`high_field_driving_force: electric_field`. No default configuration changes.

The committed template is:

`reference_tcad/pn2d_sentaurus2018_minimal6/vela/pn2d_minimal6_sweep_element_edge_gss_laux_candidate.json`

## RED and GREEN

The first attempt was rejected before Newton because the Release
`vela_example_runner.exe` predated the committed element-edge parser support.
The source parser already accepted the mode, so no source patch was needed.
Rebuilding the Release runner converted the same configuration from RED to
GREEN.

The stale-runner RED evidence is retained at:

`build-release/pn2d-minimal6-element-edge-self-consistent-20260725-stale-runner-red`

## Deterministic sweep evidence

Candidate sweep roots:

- `build-release/pn2d-minimal6-element-edge-self-consistent-20260725-a`
- `build-release/pn2d-minimal6-element-edge-self-consistent-20260725-b`

Each root contains 40 accepted exact checkpoints and zero failed transitions.
All 40 paired state CSV files are byte-identical between A and B.

Phase F comparison roots:

- `build-release/pn2d-minimal6-element-edge-phase-f-20260725-a`
- `build-release/pn2d-minimal6-element-edge-phase-f-20260725-b`

The existing independent Phase F verifier reports:

- `status: passed`
- `deterministic_pair: true`
- `verified_output_count: 6`
- `outcome: model_difference`

The more specific config-decision outcome in this report is
`config_improves_but_misses_target`: the candidate improves QFP relative to
the production baseline through the already classified electric-field
mobility setting, and greatly improves avalanche source support, but misses
the frozen full-lattice targets.

## Dependency-chain comparison

| Quantity | Production baseline median | Electric mobility median | Element-edge median | Frozen target |
| --- | ---: | ---: | ---: | ---: |
| electrostatic potential (V) | 5.37e-12 | 5.37e-12 | 5.37e-12 | 1e-6 max |
| electron QFP (V) | 0.0548737 | 0.0514769 | 0.0514769 | 0.01 |
| hole QFP (V) | 0.0620671 | 0.0591186 | 0.0591186 | 0.01 |
| electron density (dex) | 0.921833 | 0.864770 | 0.864770 | 0.10 |
| hole density (dex) | 1.042676 | 0.993144 | 0.993144 | 0.10 |
| electron element mobility (dex) | 0.465174 | 0.465257 | 0.465257 | diagnostic |
| hole element mobility (dex) | 0.292739 | 0.292733 | 0.292733 | diagnostic |
| electron directed current (dex) | 0.823474 | 0.823152 | 0.823152 | 0.10 |
| hole directed current (dex) | 0.939922 | 0.939636 | 0.939636 | 0.10 |
| terminal current (dex) | 0.647139 | 0.646672 | 0.646672 | 0.10 |
| impact source (dex) | 10.7570 | 10.7649 | 1.23508 | 0.30 |

The element-edge source change is downstream of the self-consistent QFP
branch. It does not repair the earlier QFP, density, mobility, or terminal
current mismatch. It does isolate and remove approximately 9.53 dex of the
production triangle-source discrepancy.

## Avalanche-active bias range

The two topologies agree after the verified cell permutation. The median
impact-source errors over the two topologies are:

| Bias (V) | Impact source error (dex) |
| ---: | ---: |
| -10 | 1.39034 |
| -11 | 1.07982 |
| -12 | 0.843838 |
| -13 | 0.660265 |
| -14 | 0.514620 |
| -15 | 0.397109 |
| -16 | 0.300915 |
| -17 | 0.221172 |
| -18 | 0.154333 |
| -19 | 0.0977667 |
| -20 | 0.0494843 |

The -16 V through -20 V median is `0.154333 dex`. The full 40-state median is
`1.23508 dex` because the -1 V through -9 V native sources approach numerical
zero; their ratios are useful as low-signal diagnostics but do not represent
the avalanche-active regime.

## Conservation and nonlinear behavior

At the 40 exact checkpoints:

- Newton iterations: minimum 1, median 1, maximum 1
- terminal SG versus residual current maximum relative difference:
  `3.36e-16`
- anode plus cathode KCL median absolute error:
  `2.08e-24 A/um`
- anode plus cathode KCL maximum absolute error:
  `2.46e-19 A/um`
- KCL relative P95: `4.36e-6`
- KCL relative maximum: `0.01490` at the tiny-current mirror -1 V state

The candidate therefore introduces no continuation failure, terminal-current
method inconsistency, or high-bias KCL regression.

## Decision

Keep the operator opt-in. Do not change the production default in Task 9.

The evidence supports the element-edge GSS/Laux support and source mapping as
the correct Minimal6 avalanche-active diagnostic. It does not establish
general-mesh production readiness, and it does not solve the upstream
self-consistent QFP/current model difference. Task 10 must complete
general-mesh tests, full Release regression, independent scientific and code
review, and the final default decision.
